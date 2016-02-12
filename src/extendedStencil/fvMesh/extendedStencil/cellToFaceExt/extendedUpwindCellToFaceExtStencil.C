/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
 2015-11-17 AtmosFOAM, Hilary Weller, University of Reading added support for
 cyclic boundaries
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "mapDistribute.H"
#include "extendedUpwindCellToFaceExtStencil.H"
#include "cellToFaceExtStencil.H"
#include "globalMeshData.H"
#include "syncTools.H"
#include "SortableList.H"
#include "dummyTransform.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::extendedUpwindCellToFaceExtStencil::selectOppositeFaces
(
    const boolList& nonEmptyFace,
    const scalar minOpposedness,
    const label faceI,
    const label cellI,
    DynamicList<label>& oppositeFaces
) const
{
    const vectorField& areas = mesh_.faceAreas();
    const labelList& own = mesh_.faceOwner();
    const cell& cFaces = mesh_.cells()[cellI];

    SortableList<scalar> opposedness(cFaces.size(), -GREAT);

    // Pick up all the faces that oppose this one.
    forAll(cFaces, i)
    {
        label otherFaceI = cFaces[i];

        if (otherFaceI != faceI && nonEmptyFace[otherFaceI])
        {
            if ((own[otherFaceI] == cellI) == (own[faceI] == cellI))
            {
                opposedness[i] = -(areas[otherFaceI] & areas[faceI]);
            }
            else
            {
                opposedness[i] = (areas[otherFaceI] & areas[faceI]);
            }
        }
    }

    label sz = opposedness.size();

    oppositeFaces.clear();

    scalar myAreaSqr = magSqr(areas[faceI]);

    if (myAreaSqr > VSMALL)
    {
        forAll(opposedness, i)
        {
            opposedness[i] /= myAreaSqr;
        }
        // Sort in incrementing order
        opposedness.sort();

        // Pick largest no matter what
        oppositeFaces.append(cFaces[opposedness.indices()[sz-1]]);

        for (label i = sz-2; i >= 0; --i)
        {
            if (opposedness[i] < minOpposedness)
            {
                break;
            }
            oppositeFaces.append(cFaces[opposedness.indices()[i]]);
        }
    }
    else
    {
        // Sort in incrementing order
        opposedness.sort();

        // Tiny face. Do what?
        // Pick largest no matter what
        oppositeFaces.append(cFaces[opposedness.indices()[sz-1]]);
    }
}


void Foam::extendedUpwindCellToFaceExtStencil::transportStencil
(
    const boolList& nonEmptyFace,
    const labelListList& faceStencil,
    const labelListList& transformedFaceStencil,
    const scalar minOpposedness,
    const label faceI,
    const label cellI,
    const bool stencilHasNeighbour,

    DynamicList<label>& oppositeFaces,
    labelHashSet& faceStencilSet,
    labelHashSet& transformedFaceStencilSet,
    labelList& transportedStencil,
    labelList& transportedTransformedStencil
) const
{
    label globalOwn = faceStencil[faceI][0];
    label globalNei = -1;
    if (stencilHasNeighbour && faceStencil[faceI].size() >= 2)
    {
        globalNei = faceStencil[faceI][1];
    }

    selectOppositeFaces
    (
        nonEmptyFace,
        minOpposedness,
        faceI,
        cellI,
        oppositeFaces
    );

    // Collect all stencils of opposite faces
    faceStencilSet.clear();
    forAll(oppositeFaces, i)
    {
        const labelList& fStencil = faceStencil[oppositeFaces[i]];

        forAll(fStencil, j)
        {
            label globalI = fStencil[j];

            if (globalI != globalOwn && globalI != globalNei)
            {
                faceStencilSet.insert(globalI);
            }
        }
        
        const labelList& fTransformedStencil
             = transformedFaceStencil[oppositeFaces[i]];
        
        forAll(fTransformedStencil, j)
        {
            transformedFaceStencilSet.insert(fTransformedStencil[j]);
        }
    }

    // Add my owner and neighbour first.
    if (stencilHasNeighbour)
    {
        transportedStencil.setSize(faceStencilSet.size()+2);
        transportedTransformedStencil.setSize(transformedFaceStencil.size());
        label n = 0;
        transportedStencil[n++] = globalOwn;
        transportedStencil[n++] = globalNei;

        forAllConstIter(labelHashSet, faceStencilSet, iter)
        {
            if (iter.key() != globalOwn && iter.key() != globalNei)
            {
                transportedStencil[n++] = iter.key();
            }
        }
        if (n != transportedStencil.size())
        {
            FatalErrorIn
            (
                "extendedUpwindCellToFaceExtStencil::transportStencil(..)"
            )   << "problem:" << faceStencilSet
                << abort(FatalError);
        }
        n = 0;
        forAllConstIter(labelHashSet, transformedFaceStencilSet, iter)
        {
            transportedTransformedStencil[n++] = iter.key();
        }
        if (n != transportedTransformedStencil.size())
        {
            FatalErrorIn
            (
                "extendedUpwindCellToFaceExtStencil::transportStencil(..)"
            )   << "problem:" << transformedFaceStencilSet
                << abort(FatalError);
        }
    }
    else
    {
        transportedStencil.setSize(faceStencilSet.size()+1);
        transportedTransformedStencil.setSize(transformedFaceStencil.size());
        label n = 0;
        transportedStencil[n++] = globalOwn;

        forAllConstIter(labelHashSet, faceStencilSet, iter)
        {
            if (iter.key() != globalOwn)
            {
                transportedStencil[n++] = iter.key();
            }
        }
        if (n != transportedStencil.size())
        {
            FatalErrorIn
            (
                "extendedUpwindCellToFaceStencil::transportStencil(..)"
            )   << "problem:" << faceStencilSet
                << abort(FatalError);
        }
        n = 0;
        forAllConstIter(labelHashSet, transformedFaceStencilSet, iter)
        {
            transportedTransformedStencil[n++] = iter.key();
        }
        if (n != transportedTransformedStencil.size())
        {
            FatalErrorIn
            (
                "extendedUpwindCellToFaceExtStencil::transportStencil(..)"
            )   << "problem:" << transformedFaceStencilSet
                << abort(FatalError);
        }
    }
}


void Foam::extendedUpwindCellToFaceExtStencil::transportStencils
(
    const labelListList& faceStencil,
    const List<labelPairList>& transformedStencil,
    const scalar minOpposedness,
    labelListList& ownStencil,
    labelListList& neiStencil,
    List<labelPairList>& ownTransformedElements,
    List<labelPairList>& neiTransformedElements
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const label nBnd = mesh_.nFaces()-mesh_.nInternalFaces();
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    // Work arrays
    DynamicList<label> oppositeFaces;
    labelHashSet faceStencilSet;
    labelHashSet transformedFaceStencilSet;

    // For quick detection of empty faces
    boolList nonEmptyFace(mesh_.nFaces(), true);
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (isA<emptyPolyPatch>(pp))
        {
            label faceI = pp.start();
            forAll(pp, i)
            {
                nonEmptyFace[faceI++] = false;
            }
        }
    }


    // Do the owner side
    // ~~~~~~~~~~~~~~~~~
    // stencil is synchronised at entry so no need to swap.

    ownStencil.setSize(mesh_.nFaces());
    ownTransformedElements.setSize(mesh_.nFaces());

    // Internal faces
    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        // Get stencil as owner + neighbour + stencil from 'opposite' faces
        transportStencil
        (
            nonEmptyFace,
            faceStencil,
            transformedStencil,
            minOpposedness,
            faceI,
            own[faceI],
            true,                   //stencilHasNeighbour
            oppositeFaces,
            faceStencilSet,
            transformedFaceStencilSet,
            ownStencil[faceI],
            ownTransformedElements[faceI]
        );
    }
    // Boundary faces
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        label faceI = pp.start();

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                transportStencil
                (
                    nonEmptyFace,
                    faceStencil,
                    transformedStencil,
                    minOpposedness,
                    faceI,
                    own[faceI],
                    true,                   //stencilHasNeighbour

                    oppositeFaces,
                    faceStencilSet,
                    transformedFaceStencilSet,
                    ownStencil[faceI],
                    ownTransformedElements[faceI]
                );
                faceI++;
            }
        }
        else if (!isA<emptyPolyPatch>(pp))
        {
            forAll(pp, i)
            {
                // faceStencil does not contain neighbour
                transportStencil
                (
                    nonEmptyFace,
                    faceStencil,
                    transformedStencil,
                    minOpposedness,
                    faceI,
                    own[faceI],
                    false,                  //stencilHasNeighbour

                    oppositeFaces,
                    faceStencilSet,
                    transformedFaceStencilSet,
                    ownStencil[faceI],
                    ownTransformedElements[faceI]
                );
                faceI++;
            }
        }
    }


    // Swap coupled boundary stencil
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // No idea how to do this for the coupled stencils

    labelListList neiBndStencil(nBnd);
    for (label faceI = mesh_.nInternalFaces(); faceI < mesh_.nFaces(); faceI++)
    {
        neiBndStencil[faceI-mesh_.nInternalFaces()] = ownStencil[faceI];
    }
    //syncTools::swapBoundaryFaceList(mesh_, neiBndStencil);
    syncTools::syncBoundaryFaceList
    (
        mesh_,
        neiBndStencil,
        eqOp<labelList>(),
        dummyTransform()
    );



    // Do the neighbour side
    // ~~~~~~~~~~~~~~~~~~~~~
    // - internal faces : get opposite faces on neighbour side
    // - boundary faces : empty
    // - coupled faces  : in neiBndStencil

    neiStencil.setSize(mesh_.nFaces());
    neiTransformedElements.setSize(mesh_.nFaces());

    // Internal faces
    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        transportStencil
        (
            nonEmptyFace,
            faceStencil,
            transformedStencil,
            minOpposedness,
            faceI,
            nei[faceI],
            true,                   //stencilHasNeighbour

            oppositeFaces,
            faceStencilSet,
            transformedFaceStencilSet,
            neiStencil[faceI],
            neiTransformedElements[faceI]
        );
    }

    // Boundary faces
    // No idea how to do this for the coupled stencils
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        label faceI = pp.start();

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                neiStencil[faceI].transfer
                (
                    neiBndStencil[faceI-mesh_.nInternalFaces()]
                );
                faceI++;
            }
        }
        else
        {
            // Boundary has empty neighbour stencil
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extendedUpwindCellToFaceExtStencil::extendedUpwindCellToFaceExtStencil
(
    const cellToFaceExtStencil& stencil,
    const bool pureUpwind,
    const scalar minOpposedness
)
:
    extendedCellToFaceExtStencil(stencil.mesh()),
    pureUpwind_(pureUpwind)
{
    // Transport centred stencil to upwind/downwind face
    // upwind part stored in ownStencil_, ownTransElems,
    // downwind part stored in neiStencil_, neiTransElems

    List<labelPairList> ownTransElems;
    List<labelPairList> neiTransElems;
    transportStencils
    (
        stencil.untransformedElements(),
        stencil.transformedElements(),
        minOpposedness,
        ownStencil_,
        neiStencil_,
        ownTransElems,
        neiTransElems
    );

    // Convert owner side untransformed and transformed cell indices into
    // a schedule
    {
        List<Map<label> > compactMap(Pstream::nProcs());
        ownMapPtr_.reset
        (
            new mapDistribute
            (
                stencil.globalNumbering(),
                ownStencil_,
                stencil.mesh().globalData().globalTransforms(),
                ownTransElems,
                ownTransformedElements_,
                compactMap
            )
        );
    }

    // Convert neighbour side untransformed and transformed cell indices into
    // a schedule
    {
        List<Map<label> > compactMap(Pstream::nProcs());
        neiMapPtr_.reset
        (
            new mapDistribute
            (
                stencil.globalNumbering(),
                neiStencil_,
                stencil.mesh().globalData().globalTransforms(),
                neiTransElems,
                neiTransformedElements_,
                compactMap
            )
        );
    }

//TDB
//    // stencil now in compact form
//    // Haven't attempted this yet for transformed stencils
//    if (pureUpwind_)
//    {
//        const fvMesh& mesh = dynamic_cast<const fvMesh&>(stencil.mesh());
//
//        List<List<point> > stencilPoints(ownStencil_.size());
//
//        // Owner stencil
//        // ~~~~~~~~~~~~~
//
//        collectData(ownMapPtr_(), ownStencil_, mesh.C(), stencilPoints);
//
//        // Mask off all stencil points on wrong side of face
//        forAll(stencilPoints, faceI)
//        {
//            const point& fc = mesh.faceCentres()[faceI];
//            const vector& fArea = mesh.faceAreas()[faceI];
//
//            const List<point>& points = stencilPoints[faceI];
//            const labelList& stencil = ownStencil_[faceI];
//
//            DynamicList<label> newStencil(stencil.size());
//            forAll(points, i)
//            {
//                if (((points[i]-fc) & fArea) < 0)
//                {
//                    newStencil.append(stencil[i]);
//                }
//            }
//            if (newStencil.size() != stencil.size())
//            {
//                ownStencil_[faceI].transfer(newStencil);
//            }
//        }
//
//
//        // Neighbour stencil
//        // ~~~~~~~~~~~~~~~~~
//
//        collectData(neiMapPtr_(), neiStencil_, mesh.C(), stencilPoints);
//
//        // Mask off all stencil points on wrong side of face
//        forAll(stencilPoints, faceI)
//        {
//            const point& fc = mesh.faceCentres()[faceI];
//            const vector& fArea = mesh.faceAreas()[faceI];
//
//            const List<point>& points = stencilPoints[faceI];
//            const labelList& stencil = neiStencil_[faceI];
//
//            DynamicList<label> newStencil(stencil.size());
//            forAll(points, i)
//            {
//                if (((points[i]-fc) & fArea) > 0)
//                {
//                    newStencil.append(stencil[i]);
//                }
//            }
//            if (newStencil.size() != stencil.size())
//            {
//                neiStencil_[faceI].transfer(newStencil);
//            }
//        }
//
//        // Note: could compact schedule as well. for if cells are not needed
//        // across any boundary anymore. However relatively rare.
//    }
//}

// No idea how to modify this constructor for tranformed stencils
Foam::extendedUpwindCellToFaceExtStencil::extendedUpwindCellToFaceExtStencil
(
    const cellToFaceExtStencil& stencil
)
:
    extendedCellToFaceExtStencil(stencil.mesh()),
    pureUpwind_(true)
{
    const fvMesh& mesh = dynamic_cast<const fvMesh&>(stencil.mesh());

    // Calculate stencil points with full stencil
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    List<List<point> > stencilPoints(ownStencil_.size());
    {
        labelListList untransfoIndices(stencil.untransformedElements());
        labelListList trafoIndices;
        List<Map<label> > compactMap(Pstream::nProcs());
        mapDistribute map
        (
            stencil.globalNumbering(),
            untransfoIndices,
            stencil.mesh().globalData().globalTransforms(),
            stencil.transformedElements(),
            trafoIndices,
            compactMap
        );
        collectData
        (
            map,
            untransfoIndices,
            trafoIndices,
            mesh.C(),
            stencilPoints,
            mapDistribute::transformPosition()
        );
    }


    // Split stencil into owner and neighbour
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ownStencil_.setSize(mesh.nFaces());
    neiStencil_.setSize(mesh.nFaces());

    forAll(stencilPoints, faceI)
    {
        const point& fc = mesh.faceCentres()[faceI];
        const vector& fArea = mesh.faceAreas()[faceI];

        const List<point>& points = stencilPoints[faceI];

        const labelList& untrafo =
            stencil.untransformedElements()[faceI];
        const labelPairList& trafo =
            stencil.transformedElements()[faceI];

        DynamicList<label> newOwnStencil(untrafo.size());
        DynamicList<labelPair> newOwnTrafoStencil(
        DynamicList<label> newNeiStencil(untrafo.size());
        forAll(points, i)
        {
            if (((points[i]-fc) & fArea) > 0)
            {
                newNeiStencil.append(stencil[i]);
            }
            else
            {
                newOwnStencil.append(stencil[i]);
            }
        }
        if (newNeiStencil.size() > 0)
        {
            ownStencil_[faceI].transfer(newOwnStencil);
            neiStencil_[faceI].transfer(newNeiStencil);
        }
    }

    // Should compact schedule. Or have both return the same schedule.
    neiMapPtr_.reset(new mapDistribute(ownMapPtr_()));
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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
    const List<labelPairList>& transformedFaceStencil,
    const scalar minOpposedness,
    const label faceI,
    const label cellI,
    const cellToFaceExtStencil::neighbourLocation nbrStat,

    // Work storage
    DynamicList<label>& oppositeFaces,
    labelHashSet& faceStencilSet,
    HashSet<labelPair, labelPair::Hash<>>& transformedFaceStencilSet,
    labelList& transportedStencil,
    labelPairList& transportedTransformedStencil
) const
{
    const globalIndexAndTransform& globalTransforms =
        mesh_.globalData().globalTransforms();
    const label nullIndex = globalTransforms.nullTransformIndex();

    // Extract owner (always in untransformed part) and neighbour
    // (depends on nbrStat)
    label globalOwn = faceStencil[faceI][0];
    label globalNei = -1;
    labelPair globalTrafoNei
    (
        globalIndexAndTransform::encode
        (
            -1,
            -1,
            nullIndex
        )
    );

    if (nbrStat == cellToFaceExtStencil::UNTRANSFORMED)
    {
        globalNei = faceStencil[faceI][1];
    }
    else if (nbrStat == cellToFaceExtStencil::TRANSFORMED)
    {
        globalTrafoNei = transformedFaceStencil[faceI][0];
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
    transformedFaceStencilSet.clear();
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
        
        const labelPairList& fTransformedStencil
             = transformedFaceStencil[oppositeFaces[i]];
        
        forAll(fTransformedStencil, j)
        {
            const labelPair& globalI = fTransformedStencil[j];
            if (globalI != globalTrafoNei)
            {
                transformedFaceStencilSet.insert(globalI);
            }
        }
    }

    // Add my owner and neighbour first and then rest of stencils

    label untrafoI = 0;
    label trafoI = 0;

    if (nbrStat == cellToFaceExtStencil::NONE)
    {
        transportedStencil.setSize(faceStencilSet.size()+1);
        transportedTransformedStencil.setSize(transformedFaceStencil.size());
        transportedStencil[untrafoI++] = globalOwn;
    }
    else if (nbrStat == cellToFaceExtStencil::UNTRANSFORMED)
    {
        transportedStencil.setSize(faceStencilSet.size()+2);
        transportedTransformedStencil.setSize(transformedFaceStencil.size());
        transportedStencil[untrafoI++] = globalOwn;
        transportedStencil[untrafoI++] = globalNei;
    }
    else    // TRANSFORMED
    {
        transportedStencil.setSize(faceStencilSet.size()+1);
        transportedTransformedStencil.setSize(transformedFaceStencil.size()+1);

        transportedStencil[untrafoI++] = globalOwn;
        transportedTransformedStencil[trafoI++] = globalTrafoNei;
    }

    forAllConstIter(labelHashSet, faceStencilSet, iter)
    {
        transportedStencil[untrafoI++] = iter.key();
    }

    typedef HashSet<labelPair, labelPair::Hash<> > labelPairHashSet;

    forAllConstIter
    (
        labelPairHashSet,
        transformedFaceStencilSet,
        iter
    )
    {
        transportedTransformedStencil[trafoI++] = iter.key();
    }
}


void Foam::extendedUpwindCellToFaceExtStencil::transportStencils
(
    const labelListList& faceStencil,
    const List<labelPairList>& transformedStencil,
    const scalar minOpposedness,
    labelListList& ownStencil,
    List<labelPairList>& ownTransformedElements,
    labelListList& neiStencil,
    List<labelPairList>& neiTransformedElements
) const
{
    const globalIndexAndTransform& globalTransforms =
        mesh_.globalData().globalTransforms();
    const label nullIndex = globalTransforms.nullTransformIndex();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const label nBnd = mesh_.nFaces()-mesh_.nInternalFaces();
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    // Work arrays
    DynamicList<label> oppositeFaces;
    labelHashSet faceStencilSet;
    HashSet<labelPair, labelPair::Hash<>> transformedFaceStencilSet;

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
            cellToFaceExtStencil::UNTRANSFORMED,   // nbr in untransformed part
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
            const labelPair& transSign =
                globalTransforms.patchTransformSign()[patchI];

            if (transSign.first() == nullIndex)
            {
                // Colocated coupled patch. Do exactly as if internal face
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
                        cellToFaceExtStencil::UNTRANSFORMED,

                        oppositeFaces,
                        faceStencilSet,
                        transformedFaceStencilSet,
                        ownStencil[faceI],
                        ownTransformedElements[faceI]
                    );
                    faceI++;
                }
            }
            else
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
                        cellToFaceExtStencil::TRANSFORMED, // nbr in transfrmd part

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
                    cellToFaceExtStencil::NONE,

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
    List<labelPairList> neiBndTrafoStencil(nBnd);
    for (label faceI = mesh_.nInternalFaces(); faceI < mesh_.nFaces(); faceI++)
    {
        label bFaceI = faceI-mesh_.nInternalFaces();
        neiBndStencil[bFaceI] = ownStencil[faceI];
        neiBndTrafoStencil[bFaceI] = ownTransformedElements[faceI];
    }
    syncTools::syncBoundaryFaceList
    (
        mesh_,
        neiBndStencil,
        eqOp<labelList>(),
        dummyTransform()
    );
    syncTools::syncBoundaryFaceList
    (
        mesh_,
        neiBndTrafoStencil,
        eqOp<labelPairList>(),
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
            cellToFaceExtStencil::UNTRANSFORMED,   // nbr in untransformed part

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
                label bFaceI = faceI-mesh_.nInternalFaces();
                neiStencil[faceI].transfer(neiBndStencil[bFaceI]);
                neiTransformedElements[faceI].transfer
                (
                    neiTransformedElements[bFaceI]
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
    // upwind part stored in ownUntransformedElements_, ownTransformedElements_,
    // downwind part stored in
    // neiUntransformedElements_, neiTransformedElements_

    List<labelPairList> ownTransElems;
    List<labelPairList> neiTransElems;
    transportStencils
    (
        stencil.untransformedElements(),
        stencil.transformedElements(),
        minOpposedness,
        ownUntransformedElements_,
        ownTransElems,
        neiUntransformedElements_,
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
                ownUntransformedElements_,
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
                neiUntransformedElements_,
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
//        List<List<point> > stencilPoints(ownUntransformedElements_.size());
//
//        // Owner stencil
//        // ~~~~~~~~~~~~~
//
//        collectData(ownMapPtr_(), ownUntransformedElements_, mesh.C(), stencilPoints);
//
//        // Mask off all stencil points on wrong side of face
//        forAll(stencilPoints, faceI)
//        {
//            const point& fc = mesh.faceCentres()[faceI];
//            const vector& fArea = mesh.faceAreas()[faceI];
//
//            const List<point>& points = stencilPoints[faceI];
//            const labelList& stencil = ownUntransformedElements_[faceI];
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
//                ownUntransformedElements_[faceI].transfer(newStencil);
//            }
//        }
//
//
//        // Neighbour stencil
//        // ~~~~~~~~~~~~~~~~~
//
//        collectData(neiMapPtr_(), neiUntransformedElements_, mesh.C(), stencilPoints);
//
//        // Mask off all stencil points on wrong side of face
//        forAll(stencilPoints, faceI)
//        {
//            const point& fc = mesh.faceCentres()[faceI];
//            const vector& fArea = mesh.faceAreas()[faceI];
//
//            const List<point>& points = stencilPoints[faceI];
//            const labelList& stencil = neiUntransformedElements_[faceI];
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
//                neiUntransformedElements_[faceI].transfer(newStencil);
//            }
//        }
//
//        // Note: could compact schedule as well. for if cells are not needed
//        // across any boundary anymore. However relatively rare.
//    }
}


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

    List<List<point> > stencilPoints(ownUntransformedElements_.size());
    {
        labelListList untransfoIndices(stencil.untransformedElements());
        labelListList trafoIndices;
        List<Map<label> > compactMap(Pstream::nProcs());
        ownMapPtr_.reset
        (
            new mapDistribute
            (
                stencil.globalNumbering(),
                untransfoIndices,
                stencil.mesh().globalData().globalTransforms(),
                stencil.transformedElements(),
                trafoIndices,
                compactMap
            )
        );
        collectData
        (
            ownMapPtr_(),
            untransfoIndices,
            trafoIndices,
            mesh.C(),
            stencilPoints,
            mapDistribute::transformPosition()
        );
    }


    // Split stencil into owner and neighbour
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ownUntransformedElements_.setSize(mesh.nFaces());
    List<labelPairList> ownTransElems(mesh.nFaces());

    neiUntransformedElements_.setSize(mesh.nFaces());
    List<labelPairList> neiTransElems(mesh.nFaces());

    forAll(stencilPoints, faceI)
    {
        const point& fc = mesh.faceCentres()[faceI];
        const vector& fArea = mesh.faceAreas()[faceI];

        const List<point>& points = stencilPoints[faceI];

        const labelList& untrafo = stencil.untransformedElements()[faceI];
        const labelPairList& trafo = stencil.transformedElements()[faceI];

        DynamicList<label> newOwnStencil(untrafo.size());
        DynamicList<labelPair> newOwnTrafoStencil(trafo.size());

        DynamicList<label> newNeiStencil(untrafo.size());
        DynamicList<labelPair> newNeiTrafoStencil(trafo.size());

        // Untransformed first
        label compactI = 0;
        forAll(untrafo, i)
        {
            const point& pt = points[compactI++];

            if (((pt-fc) & fArea) > 0)
            {
                newNeiStencil.append(untrafo[i]);
            }
            else
            {
                newOwnStencil.append(untrafo[i]);
            }
        }
        // Transformed second
        forAll(trafo, i)
        {
            const point& pt = points[compactI++];

            if (((pt-fc) & fArea) > 0)
            {
                newNeiTrafoStencil.append(trafo[i]);
            }
            else
            {
                newOwnTrafoStencil.append(trafo[i]);
            }
        }


        ownUntransformedElements_[faceI].transfer(newOwnStencil);
        ownTransElems[faceI].transfer(newOwnTrafoStencil);

        neiUntransformedElements_[faceI].transfer(newNeiStencil);
        neiTransElems[faceI].transfer(newNeiTrafoStencil);
    }


    // Convert owner side untransformed and transformed cell indices into
    // a schedule
    {
        List<Map<label> > compactMap(Pstream::nProcs());
        ownMapPtr_.reset
        (
            new mapDistribute
            (
                stencil.globalNumbering(),
                ownUntransformedElements_,
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
                neiUntransformedElements_,
                stencil.mesh().globalData().globalTransforms(),
                neiTransElems,
                neiTransformedElements_,
                compactMap
            )
        );
    }
}


// ************************************************************************* //

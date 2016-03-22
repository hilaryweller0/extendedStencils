/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "CFCCellToCellExtStencil.H"
#include "syncTools.H"
#include "SortableList.H"
#include "emptyPolyPatch.H"
#include "cellToFaceExtStencil.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::CFCCellToCellExtStencil::calcFaceBoundaryData
(
    labelList& neiGlobal,
    labelPairList& neiTrafoGlobal
) const
{
    const globalIndexAndTransform& globalTransforms =
        mesh().globalData().globalTransforms();
    const label nullIndex = globalTransforms.nullTransformIndex();
    const polyBoundaryMesh& patches = mesh().boundaryMesh();
    const label nBnd = mesh().nFaces()-mesh().nInternalFaces();

    neiGlobal.setSize(nBnd);
    neiGlobal = -1;
    neiTrafoGlobal.setSize(nBnd);
    neiTrafoGlobal = labelPair(-1, -1);

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            const labelPair& transSign =
                globalTransforms.patchTransformSign()[patchI];

            const labelUList& fc = pp.faceCells();

            labelList globalCells(fc.size());
            forAll(fc, i)
            {
                globalCells[i] = globalNumbering().toGlobal(fc[i]);
            }

            if (transSign.first() == nullIndex)
            {
                // Get owner cells (untransformed)
                forAll(pp, i)
                {
                    label bFaceI = pp.start()+i-mesh().nInternalFaces();
                    neiGlobal[bFaceI] = globalCells[i];
                }
            }
            else
            {
                // Get owner cells (transformed)
                labelPairList trafoCells
                (
                    cellToFaceExtStencil::transform
                    (
                        mesh(),
                        patchI,
                        globalNumbering(),
                        globalCells,
                        true        // Sending patch
                    )
                );

                forAll(pp, i)
                {
                    label bFaceI = pp.start()+i-mesh().nInternalFaces();
                    neiTrafoGlobal[bFaceI] = trafoCells[i];
                }
            }
        }
        else if (!isA<emptyPolyPatch>(pp))
        {
            // For noncoupled faces get the boundary face.
            forAll(pp, i)
            {
                label bFaceI = pp.start()+i-mesh().nInternalFaces();
                neiGlobal[bFaceI] =
                    globalNumbering().toGlobal(mesh().nCells()+bFaceI);
            }
        }
    }
    syncTools::swapBoundaryFaceList(mesh(), neiGlobal);

    // Note: cannot do
    //    syncTools::swapBoundaryFaceList(mesh(), neiTrafoGlobal);
    // since uses SubField and Field ops not definable for a Pair so workaround
    {
        labelList a(neiTrafoGlobal.size());
        labelList b(neiTrafoGlobal.size());
        forAll(neiTrafoGlobal, i)
        {
            a[i] = neiTrafoGlobal[i][0];
            b[i] = neiTrafoGlobal[i][1];
        }
        syncTools::swapBoundaryFaceList(mesh(), a);
        syncTools::swapBoundaryFaceList(mesh(), b);
        forAll(neiTrafoGlobal, i)
        {
            neiTrafoGlobal[i][0] = a[i];
            neiTrafoGlobal[i][1] = b[i];
        }
    }
}


void Foam::CFCCellToCellExtStencil::calcCellStencil
(
    labelListList& globalCellCells,
    List<labelPairList>& globalTrafoCellCells
) const
{
    const label nBnd = mesh().nFaces()-mesh().nInternalFaces();
    const labelList& own = mesh().faceOwner();
    const labelList& nei = mesh().faceNeighbour();


    // Calculate coupled neighbour (in global numbering)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList neiGlobal(nBnd);
    labelPairList neiTrafoGlobal(nBnd);
    calcFaceBoundaryData(neiGlobal, neiTrafoGlobal);


    // Determine cellCells in global numbering
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    globalCellCells.setSize(mesh().nCells());
    globalTrafoCellCells.setSize(mesh().nCells());

    forAll(globalCellCells, cellI)
    {
        const cell& cFaces = mesh().cells()[cellI];

        label nTrafo = 0;
        label nUntrafo = 0;
        forAll(cFaces, i)
        {
            label faceI = cFaces[i];
            label bFaceI = faceI - mesh().nInternalFaces();

            if (bFaceI >= 0)
            {
                if (neiTrafoGlobal[bFaceI] != labelPair(-1, -1))
                {
                    nTrafo++;
                }
                else if (neiGlobal[bFaceI] != -1)
                {
                    nUntrafo++;
                }
            }
            else
            {
                nUntrafo++;
            }
        }

        labelList& cCells = globalCellCells[cellI];
        cCells.setSize(nUntrafo+1);
        nUntrafo = 0;

        labelPairList& cTrafoCells = globalTrafoCellCells[cellI];
        cTrafoCells.setSize(nTrafo+1);
        nTrafo = 0;

        // Myself
        cCells[nUntrafo++] = globalNumbering().toGlobal(cellI);

        // Collect neighbouring cells/faces
        forAll(cFaces, i)
        {
            label faceI = cFaces[i];

            if (mesh().isInternalFace(faceI))
            {
                label nbrCellI = own[faceI];
                if (nbrCellI == cellI)
                {
                    nbrCellI = nei[faceI];
                }
                cCells[nUntrafo++] = globalNumbering().toGlobal(nbrCellI);
            }
            else
            {
                label bFaceI = faceI - mesh().nInternalFaces();

                if (neiTrafoGlobal[bFaceI] != labelPair(-1, -1))
                {
                    cTrafoCells[nTrafo++] = neiTrafoGlobal[bFaceI];
                }
                else if (neiGlobal[bFaceI] != -1)
                {
                    cCells[nUntrafo++] = neiGlobal[bFaceI];
                }
            }
        }
        cCells.setSize(nUntrafo);
        cTrafoCells.setSize(nTrafo);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CFCCellToCellExtStencil::CFCCellToCellExtStencil(const polyMesh& mesh)
:
    cellToCellExtStencil(mesh)
{
    // Calculate per cell the (face) connected cells (in global numbering)
    calcCellStencil(untransformedElements_, transformedElements_);
}


// ************************************************************************* //

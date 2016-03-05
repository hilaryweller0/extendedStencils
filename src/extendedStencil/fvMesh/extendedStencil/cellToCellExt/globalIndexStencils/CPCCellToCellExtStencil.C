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

#include "CPCCellToCellExtStencil.H"
#include "syncTools.H"
#include "dummyTransform.H"
#include "cellToFaceExtStencil.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Calculates per point on coupled patch the pointCells
// (with optional transformation added)
void Foam::CPCCellToCellExtStencil::calcPointBoundaryData
(
    const boolList& isValidBFace,
    const labelList& boundaryPoints,
    Map<labelList>& neiGlobal,
    Map<labelPairList>& neiTrafoGlobal
) const
{
    const globalIndexAndTransform& globalTransforms =
        mesh().globalData().globalTransforms();
    const label nullIndex = globalTransforms.nullTransformIndex();
    const polyBoundaryMesh& patches = mesh().boundaryMesh();


    // Count number of points

    label nUntrafo = 0;
    label nTrafo = 0;

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            const labelPair& transSign =
                globalTransforms.patchTransformSign()[patchI];
            if (transSign.first() == nullIndex)
            {
                nUntrafo += pp.nPoints();
            }
            else
            {
                nTrafo += pp.nPoints();
            }
        }
    }

    // Size. Assume 4 cells per point
    neiGlobal.resize(2*4*nUntrafo);
    neiTrafoGlobal.resize(2*4*nTrafo);

    labelHashSet pointGlobals;

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            const labelList& meshPoints = pp.meshPoints();

            const labelPair& transSign =
                globalTransforms.patchTransformSign()[patchI];

            forAll(meshPoints, i)
            {
                label pointI = meshPoints[i];
                // Calculate all local cells connected to the point
                labelList pCells
                (
                    calcFaceCells
                    (
                        isValidBFace,
                        mesh().pointFaces()[pointI],
                        pointGlobals
                    )
                );

                if (transSign.first() == nullIndex)
                {
                    neiGlobal.insert(pointI, pCells);
                }
                else
                {
                    // Add sending transform
                    labelPairList pTrafoCells
                    (
                        cellToFaceExtStencil::transform
                        (
                            mesh(),
                            patchI,
                            globalNumbering(),
                            pCells,
                            true        // Sending patch
                        )
                    );
                    neiTrafoGlobal.insert(pointI, pTrafoCells);
                }
            }
        }
    }


    {
        HashSet<label, Hash<label>> set;
        unionEqOp<label, Hash<label>> op(set);

        syncTools::syncPointMap
        (
            mesh(),
            neiGlobal,
            op,
            Foam::dummyTransform()      // dummy transformation
        );
    }
    {
        HashSet<labelPair, labelPair::Hash<>> set;
        unionEqOp<labelPair, labelPair::Hash<>> op(set);

        syncTools::syncPointMap
        (
            mesh(),
            neiTrafoGlobal,
            op,
            Foam::dummyTransform()      // dummy transformation
        );
    }
}


// Calculates per cell the neighbour data (= cell or boundary in global
// numbering). First element is always cell itself!
void Foam::CPCCellToCellExtStencil::calcCellStencil
(
    labelListList& globalCellCells,
    List<labelPairList>& globalTrafoCellCells
) const
{
    // Calculate points on coupled patches
    const labelList& boundaryPoints =
        mesh().globalData().coupledPatch().meshPoints();


    // Mark boundary faces to be included in stencil (i.e. not coupled or empty)
    boolList isValidBFace;
    validBoundaryFaces(isValidBFace);


    globalCellCells.setSize(mesh().nCells());
    globalTrafoCellCells.setSize(mesh().nCells());


    // 1. Do local (untranformed) pointCells
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelHashSet workSet;

    for (label pointI = 0; pointI < mesh().nPoints(); pointI++)
    {
        labelList pGlobals
        (
            calcFaceCells
            (
                isValidBFace,
                mesh().pointFaces()[pointI],
                workSet
            )
        );

        const labelList& pCells = mesh().pointCells(pointI);

        forAll(pCells, j)
        {
            label cellI = pCells[j];

            merge
            (
                globalNumbering().toGlobal(cellI),
                pGlobals,
                globalCellCells[cellI]
            );
        }
    }


    // 2. Swap pointCells for coupled points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Map<labelList> neiGlobal;
    Map<labelPairList> neiTrafoGlobal;
    calcPointBoundaryData
    (
        isValidBFace,
        boundaryPoints,
        neiGlobal,
        neiTrafoGlobal
    );


    // 3. Merge coupled data into local data
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    HashSet<labelPair, labelPair::Hash<> > trafoSet;

    forAll(boundaryPoints, i)
    {
        label pointI = boundaryPoints[i];
        const labelList& pCells = mesh().pointCells(pointI);

        // Distribute neighbour data to all pointCells

        Map<labelList>::const_iterator fnd = neiGlobal.find(pointI);
        if (fnd != neiGlobal.end())
        {
            const labelList& pGlobals = fnd();

            forAll(pCells, j)
            {
                label cellI = pCells[j];
                merge
                (
                    globalNumbering().toGlobal(cellI),
                    pGlobals,
                    globalCellCells[cellI]
                );
            }
        }

        Map<labelPairList>::const_iterator fnd2 =
            neiTrafoGlobal.find(pointI);
        if (fnd2 != neiTrafoGlobal.end())
        {
            const labelPairList& pGlobals = fnd2();

            forAll(pCells, j)
            {
                label cellI = pCells[j];

                // Merge pGlobals into globalTrafoCellCells
                trafoSet.clear();
                label sz = globalTrafoCellCells[cellI].size();
                trafoSet.resize(2*(sz+pGlobals.size()));
                trafoSet.insert(globalTrafoCellCells[cellI]);

                // Check that transformed are not present in untransformed bit
                forAll(pGlobals, pI)
                {
                    const labelPair& lp = pGlobals[pI];
                    label index = globalIndexAndTransform::index(lp);
                    label procI = globalIndexAndTransform::processor(lp);
                    label globalI =
                        globalNumbering().toGlobal(procI, index);
                    if (findIndex(globalCellCells[cellI], globalI) == -1)
                    {
                        trafoSet.insert(lp);
                    }
                }

                globalTrafoCellCells[cellI] = trafoSet.toc();
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CPCCellToCellExtStencil::CPCCellToCellExtStencil(const polyMesh& mesh)
:
    cellToCellExtStencil(mesh)
{
    // Calculate per cell the (point) connected cells (in global numbering)
    calcCellStencil(untransformedElements_, transformedElements_);
}


// ************************************************************************* //

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

#include "CECCellToCellExtStencil.H"
#include "syncTools.H"
#include "dummyTransform.H"
#include "cellToFaceExtStencil.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Calculates per edge the neighbour data (= edgeCells)
void Foam::CECCellToCellExtStencil::calcEdgeBoundaryData
(
    const boolList& isValidBFace,
    EdgeMap<labelList>& neiGlobal,
    EdgeMap<labelPairList>& neiTrafoGlobal
) const
{
    const globalIndexAndTransform& globalTransforms =
        mesh().globalData().globalTransforms();
    const label nullIndex = globalTransforms.nullTransformIndex();
    const polyBoundaryMesh& patches = mesh().boundaryMesh();


    // Count number of edges

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
                nUntrafo += pp.nEdges();
            }
            else
            {
                nTrafo += pp.nEdges();
            }
        }
    }


    // Size. Assume 2 cells per edge
    neiGlobal.resize(2*2*nUntrafo);
    neiTrafoGlobal.resize(2*2*nTrafo);

    labelHashSet edgeGlobals;

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            const labelList& meshEdges = pp.meshEdges();
            const labelPair& transSign =
                globalTransforms.patchTransformSign()[patchI];

            forAll(meshEdges, i)
            {
                label edgeI = meshEdges[i];
                const edge& e = mesh().edges()[edgeI];

                // Calculate all local cells connected to the edge
                labelList pCells
                (
                    calcFaceCells
                    (
                        isValidBFace,
                        mesh().edgeFaces(edgeI),
                        edgeGlobals
                    )
                );

                if (transSign.first() == nullIndex)
                {
                    neiGlobal.insert(e, pCells);
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
                    neiTrafoGlobal.insert(e, pTrafoCells);
                }
            }
        }
    }


    {
        HashSet<label, Hash<label>> set;
        unionEqOp<label, Hash<label>> op(set);

        syncTools::syncEdgeMap
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

        syncTools::syncEdgeMap
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
void Foam::CECCellToCellExtStencil::calcCellStencil
(
    labelListList& globalCellCells,
    List<labelPairList>& globalTrafoCellCells
) const
{
    // Calculate edges on coupled patches
    const labelList boundaryEdges
    (
         mesh().globalData().coupledPatch().meshEdges
        (
            mesh().edges(),
            mesh().pointEdges()
        )
    );

    //{
    //    OFstream str(mesh().time().path()/"boundaryEdges.obj");
    //    Pout<< "DUmping boundary edges to " << str.name() << endl;
    //
    //    label vertI = 0;
    //    forAll(boundaryEdges, i)
    //    {
    //        label edgeI = boundaryEdges[i];
    //        const edge& e = mesh().edges()[edgeI];
    //        const point& p0 = mesh().points()[e[0]];
    //        const point& p1 = mesh().points()[e[1]];
    //
    //        Pout<< "boundary edge " << edgeI << " between " << p0 << p1
    //            << endl;
    //
    //        meshTools::writeOBJ(str, p0);
    //        vertI++;
    //        meshTools::writeOBJ(str, p1);
    //        vertI++;
    //        str << "l " << vertI-1 << ' ' << vertI << nl;
    //    }
    //}


    // Mark boundary faces to be included in stencil (i.e. not coupled or empty)
    boolList isValidBFace;
    validBoundaryFaces(isValidBFace);


    globalCellCells.setSize(mesh().nCells());
    globalTrafoCellCells.setSize(mesh().nCells());


    // 1. Do local (untranformed) edgeCells
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelHashSet workSet;

    for (label edgeI = 0; edgeI < mesh().nEdges(); edgeI++)
    {
        labelList eGlobals
        (
            calcFaceCells
            (
                isValidBFace,
                mesh().edgeFaces(edgeI),
                workSet
            )
        );

        const labelList& eCells = mesh().edgeCells(edgeI);

        forAll(eCells, j)
        {
            label cellI = eCells[j];

            merge
            (
                globalNumbering().toGlobal(cellI),
                eGlobals,
                globalCellCells[cellI]
            );
        }
    }



    // 2. Swap edgeCells for coupled edges
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Swap edgeCells for coupled edges. Note: use EdgeMap for now since we've
    // got syncTools::syncEdgeMap for those. Should be replaced with Map and
    // syncTools functionality to handle those.
    EdgeMap<labelList> neiGlobal;
    EdgeMap<labelPairList> neiTrafoGlobal;
    calcEdgeBoundaryData
    (
        isValidBFace,
        neiGlobal,
        neiTrafoGlobal
    );


    // 3. Merge coupled data into local data
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    HashSet<labelPair, labelPair::Hash<> > trafoSet;

    forAll(boundaryEdges, i)
    {
        label edgeI = boundaryEdges[i];
        const edge& e = mesh().edges()[edgeI];
        const labelList& eCells = mesh().edgeCells(edgeI);

        EdgeMap<labelList>::const_iterator fnd = neiGlobal.find(e);
        if (fnd != neiGlobal.end())
        {
            const labelList& eGlobals = fnd();

            forAll(eCells, j)
            {
                label cellI = eCells[j];
                merge
                (
                    globalNumbering().toGlobal(cellI),
                    eGlobals,
                    globalCellCells[cellI]
                );
            }
        }

        EdgeMap<labelPairList>::const_iterator fnd2 =
            neiTrafoGlobal.find(e);
        if (fnd2 != neiTrafoGlobal.end())
        {
            const labelPairList& eGlobals = fnd2();

            forAll(eCells, j)
            {
                label cellI = eCells[j];

                // Merge eGlobals into globalTrafoCellCells
                trafoSet.clear();
                label sz = globalTrafoCellCells[cellI].size();
                trafoSet.resize(2*(sz+eGlobals.size()));
                trafoSet.insert(globalTrafoCellCells[cellI]);

                // Check that transformed are not present in untransformed bit
                forAll(eGlobals, pI)
                {
                    const labelPair& lp = eGlobals[pI];
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

Foam::CECCellToCellExtStencil::CECCellToCellExtStencil(const polyMesh& mesh)
:
    cellToCellExtStencil(mesh)
{
    // Calculate per cell the (edge) connected cells (in global numbering)
    calcCellStencil(untransformedElements_, transformedElements_);
}


// ************************************************************************* //

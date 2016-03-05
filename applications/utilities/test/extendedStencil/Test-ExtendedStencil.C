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

Application
    testExtendedStencil

Description
    Test app for determining extended stencil.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "Time.H"
//#include "mapDistribute.H"
#include "OFstream.H"
#include "meshTools.H"
//#include "FECCellToFaceStencil.H"
//#include "CFCCellToFaceStencil.H"
//#include "CPCCellToFaceStencil.H"
//#include "CECCellToFaceStencil.H"

#include "CPCCellToFaceExtStencil.H"
#include "extendedCentredCellToFaceExtStencil.H"
#include "extendedUpwindCellToFaceExtStencil.H"

//#include "centredCFCCellToFaceStencilObject.H"
//#include "centredFECCellToFaceStencilObject.H"
//#include "centredCPCCellToFaceStencilObject.H"
//#include "centredCECCellToFaceStencilObject.H"

//#include "upwindFECCellToFaceStencilObject.H"
//#include "upwindCPCCellToFaceStencilObject.H"
//#include "upwindCECCellToFaceStencilObject.H"

//#include "upwindCFCCellToFaceStencilObject.H"
//#include "centredCFCFaceToCellStencilObject.H"

#include "centredCECCellToCellStencilObject.H"
#include "centredCFCCellToCellStencilObject.H"
#include "centredCPCCellToCellStencilObject.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void writeStencilOBJ
(
    const fileName& fName,
    const point& fc,
    const List<point>& stencilCc
)
{
    OFstream str(fName);
    label vertI = 0;

    meshTools::writeOBJ(str, fc);
    vertI++;

    forAll(stencilCc, i)
    {
        meshTools::writeOBJ(str, stencilCc[i]);
        vertI++;
        str << "l 1 " << vertI << nl;
    }
}


// Stats
void writeStencilStats(const labelListList& stencil)
{
    label sumSize = 0;
    label nSum = 0;
    label minSize = labelMax;
    label maxSize = labelMin;

    forAll(stencil, i)
    {
        const labelList& sCells = stencil[i];

        if (sCells.size() > 0)
        {
            sumSize += sCells.size();
            nSum++;
            minSize = min(minSize, sCells.size());
            maxSize = max(maxSize, sCells.size());
        }
    }
    reduce(sumSize, sumOp<label>());
    reduce(nSum, sumOp<label>());
    sumSize /= nSum;

    reduce(minSize, minOp<label>());
    reduce(maxSize, maxOp<label>());

    Info<< "Stencil size :" << nl
        << "    average : " << sumSize << nl
        << "    min     : " << minSize << nl
        << "    max     : " << maxSize << nl
        << endl;
}


// Main program:

int main(int argc, char *argv[])
{
    #include "addTimeOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"

    // Get times list
    instantList Times = runTime.times();
    #include "checkTimeOptions.H"
    runTime.setTime(Times[startTime], startTime);
    #include "createMesh.H"

    // Force calculation of extended edge addressing
    const labelListList& edgeFaces = mesh.edgeFaces();
    const labelListList& edgeCells = mesh.edgeCells();
    const labelListList& pointCells = mesh.pointCells();
    Info<< "dummy:" << edgeFaces.size() + edgeCells.size() + pointCells.size()
        << endl;


    const globalIndexAndTransform& globalTransforms =
        mesh.globalData().globalTransforms();
    const label nullIndex = globalTransforms.nullTransformIndex();
    Debug(nullIndex);


//     // Centred, semi-extended stencil (edge cells only)
//     // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//     {
//         CPCCellToCellExtStencil2 cpcStencil(mesh);
//         Info<< "untransformed:" << endl;
//         writeStencilStats(cpcStencil.untransformedElements());
// 
//         for (label cellI = 0; cellI < mesh.nCells(); cellI++)
//         {
//             const labelList& cCells = cpcStencil.untransformedElements()[cellI];
//             const labelPairList& trafoCells =
//                 cpcStencil.transformedElements()[cellI];
// 
//             Pout<< "cell:" << cellI << " at:"
//                 << mesh.cellCentres()[cellI] << endl;
//             forAll(cCells, i)
//             {
//                 label index = cCells[i];
//                 if (index < mesh.nCells())
//                 {
//                     Pout<< "    cell:" << index
//                         << " at:" << mesh.cellCentres()[index]
//                         << endl;
//                 }
//                 else
//                 {
//                     label faceI = index - mesh.nCells() + mesh.nInternalFaces();
//                     Pout<< "    face:" << faceI
//                         << " at:" << mesh.faceCentres()[faceI]
//                         << endl;
//                 }
//             }
//             forAll(trafoCells, i)
//             {
//                 const labelPair& lp = trafoCells[i];
//                 label index = globalIndexAndTransform::index(lp);
//                 label procI = globalIndexAndTransform::processor(lp);
//                 label trafoI = globalIndexAndTransform::transformIndex(lp);
//                 if (index < mesh.nCells())
//                 {
//                     Pout<< "    cell:" << index
//                         << " at:" << mesh.cellCentres()[index]
//                         << " through:" << trafoI
//                         << endl;
//                 }
//                 else
//                 {
//                     label faceI = index - mesh.nCells() + mesh.nInternalFaces();
//                     Pout<< "    face:" << faceI
//                         << " at:" << mesh.faceCentres()[faceI]
//                         << " through:" << trafoI
//                         << endl;
//                 }
//             }
//         }
//     }

    // Centred, semi-extended stencil (edge cells only)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        CPCCellToFaceExtStencil cpcStencil(mesh);

        const extendedCentredCellToFaceExtStencil addressing(cpcStencil);

        Debug(addressing.untransformedElements());
        Debug(addressing.transformedElements());

        Info<< "untransformed:" << endl;
        writeStencilStats(addressing.untransformedElements());

        // Collect stencil cell centres
        List<List<point>> stencilPoints(mesh.nFaces());
        //addressing.collectData
        //(
        //    mesh.C(),
        //    stencilPoints
        //);
        extendedCellToFaceExtStencil::collectData
        (
            addressing.map(),
            addressing.untransformedElements(),
            addressing.transformedElements(),
            mesh.C(),
            stencilPoints,
            mapDistribute::transformPosition()
        );

        forAll(stencilPoints, faceI)
        {
            Debug(mesh.faceCentres()[faceI]);
            //Debug(addressing.untransformedElements()[faceI]);
            //Debug(addressing.transformedElements()[faceI]);
            Debug(stencilPoints[faceI]);
        }


       forAll(stencilPoints, faceI)
       {
           writeStencilOBJ
           (
               runTime.path()/"faceCPCCell" + Foam::name(faceI) + ".obj",
               mesh.faceCentres()[faceI],
               stencilPoints[faceI]
           );
       }
   }

    // Upwinded, semi-extended stencil (edge cells only)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        CPCCellToFaceExtStencil cpcStencil(mesh);

        const extendedUpwindCellToFaceExtStencil addressing(cpcStencil);


        // Collect stencil cell centres
        List<List<point>> ownStencilPoints(mesh.nFaces());
        addressing.collectOwnPositions(mesh.C(), ownStencilPoints);
        List<List<point>> neiStencilPoints(mesh.nFaces());
        addressing.collectNeiPositions(mesh.C(), neiStencilPoints);

        forAll(ownStencilPoints, faceI)
        {
            Debug(mesh.faceCentres()[faceI]);
            Debug(ownStencilPoints[faceI]);
            Debug(neiStencilPoints[faceI]);
        }
   }

//XXXXXX
//    // Evaluate
//    List<List<scalar>> stencilData(faceStencils.size());
//    collectStencilData
//    (
//        distMap,
//        faceStencils,
//        vf,
//        stencilData
//    );
//    for (label faci = 0; faci < mesh.nInternalFaces(); faci++)
//    {
//        const scalarList& stData = stencilData[faceI];
//        const scalarList& stWeight = fit[faceI];
//
//        forAll(stData, i)
//        {
//            sf[faceI] += stWeight[i]*stData[i];
//        }
//    }
//    See finiteVolume/lnInclude/leastSquaresGrad.C
//XXXXXX

    Pout<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Class
    Foam::upwindCPCCellToFaceExtStencilObject

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef upwindCPCCellToFaceExtStencilObject_H
#define upwindCPCCellToFaceExtStencilObject_H

#include "extendedUpwindCellToFaceExtStencil.H"
#include "CPCCellToFaceExtStencil.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class upwindCPCCellToFaceExtStencilObject Declaration
\*---------------------------------------------------------------------------*/

class upwindCPCCellToFaceExtStencilObject
:
    public MeshObject
    <
        fvMesh,
        TopologicalMeshObject,
        upwindCPCCellToFaceExtStencilObject
    >,
    public extendedUpwindCellToFaceExtStencil
{

public:

    TypeName("upwindCPCCellToFaceExtStencil");

    // Constructors

        //- Construct from uncompacted face stencil
        explicit upwindCPCCellToFaceExtStencilObject
        (
            const fvMesh& mesh
        )
        :
            MeshObject
            <
                fvMesh,
                Foam::TopologicalMeshObject,
                upwindCPCCellToFaceExtStencilObject
            >(mesh),
            extendedUpwindCellToFaceExtStencil
            (
                CPCCellToFaceExtStencil(mesh)
            )
        {
            if (extendedCellToFaceExtStencil::debug)
            {
                Info<< "Generated off-centred stencil " << type()
                    << nl << endl;
                writeStencilStats
                (
                    Info, ownUntransformedElements(), ownTransformedElements(),
                    ownMap()
                );
                writeStencilStats
                (
                    Info, neiUntransformedElements(), neiTransformedElements(),
                    neiMap()
                );
            }
        }


    //- Destructor
    virtual ~upwindCPCCellToFaceExtStencilObject()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

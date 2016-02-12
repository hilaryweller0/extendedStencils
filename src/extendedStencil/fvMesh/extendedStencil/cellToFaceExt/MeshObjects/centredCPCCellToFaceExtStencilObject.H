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
    Foam::centredCPCCellToFaceExtStencilObject

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef centredCPCCellToFaceExtStencilObject_H
#define centredCPCCellToFaceExtStencilObject_H

#include "extendedCentredCellToFaceExtStencil.H"
#include "CPCCellToFaceExtStencil.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class centredCPCCellToFaceExtStencilObject Declaration
\*---------------------------------------------------------------------------*/

class centredCPCCellToFaceExtStencilObject
:
    public MeshObject
    <
        fvMesh,
        TopologicalMeshObject,
        centredCPCCellToFaceExtStencilObject
    >,
    public extendedCentredCellToFaceExtStencil
{

public:

    TypeName("centredCPCCellToFaceExtStencil");

    // Constructors

        //- Construct from uncompacted face stencil
        explicit centredCPCCellToFaceExtStencilObject
        (
            const fvMesh& mesh
        )
        :
            MeshObject
            <
                fvMesh,
                Foam::TopologicalMeshObject,
                centredCPCCellToFaceExtStencilObject
            >(mesh),
            extendedCentredCellToFaceExtStencil(CPCCellToFaceExtStencil(mesh))
        {
            if (extendedCellToFaceExtStencil::debug)
            {
                Info<< "Generated centred stencil " << type()
                    << nl << endl;
                writeStencilStats
                (
                    Info, untransformedElements(), transformedElements(),
                    map()
                );
            }
        }


    //- Destructor
    virtual ~centredCPCCellToFaceExtStencilObject()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

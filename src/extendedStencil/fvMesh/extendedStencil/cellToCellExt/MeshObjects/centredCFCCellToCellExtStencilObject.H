/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::centredCFCCellToCellExtStencilObject

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef centredCFCCellToCellExtStencilObject_H
#define centredCFCCellToCellExtStencilObject_H

#include "extendedCentredCellToCellExtStencil.H"
#include "CFCCellToCellExtStencil.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class centredCFCCellToCellExtStencilObject Declaration
\*---------------------------------------------------------------------------*/

class centredCFCCellToCellExtStencilObject
:
    public MeshObject
    <
        fvMesh,
        TopologicalMeshObject,
        centredCFCCellToCellExtStencilObject
    >,
    public extendedCentredCellToCellExtStencil
{

public:

    TypeName("centredCFCCellToCellExtStencil");

    // Constructors

        //- Construct from uncompacted cell stencil
        explicit centredCFCCellToCellExtStencilObject
        (
            const fvMesh& mesh
        )
        :
            MeshObject
            <
                fvMesh,
                Foam::TopologicalMeshObject,
                centredCFCCellToCellExtStencilObject
            >(mesh),
            extendedCentredCellToCellExtStencil(CFCCellToCellExtStencil(mesh))
        {}


    //- Destructor
    virtual ~centredCFCCellToCellExtStencilObject()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

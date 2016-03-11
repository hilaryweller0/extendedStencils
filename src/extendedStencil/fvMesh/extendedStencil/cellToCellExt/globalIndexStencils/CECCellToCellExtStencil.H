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

Class
    Foam::CECCellToCellExtStencil

Description

SourceFiles
    CECCellToCellExtStencil.C

\*---------------------------------------------------------------------------*/

#ifndef CECCellToCellExtStencil_H
#define CECCellToCellExtStencil_H

#include "cellToCellExtStencil.H"
#include "boolList.H"
#include "HashSet.H"
#include "Map.H"
#include "EdgeMap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class CECCellToCellExtStencil Declaration
\*---------------------------------------------------------------------------*/

class CECCellToCellExtStencil
:
    public cellToCellExtStencil
{
    // Private Member Functions

        //- Calculates per edge the neighbour data (= edgeCells)
        void calcEdgeBoundaryData
        (
            const boolList& isValidBFace,
            EdgeMap<labelList>& neiGlobal,
            EdgeMap<labelPairList>& neiTrafoGlobal
        ) const;

        void calcCellStencil(labelListList&, List<labelPairList>&) const;


public:

    // Constructors

        //- Construct from all cells and boundary faces
        explicit CECCellToCellExtStencil(const polyMesh&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

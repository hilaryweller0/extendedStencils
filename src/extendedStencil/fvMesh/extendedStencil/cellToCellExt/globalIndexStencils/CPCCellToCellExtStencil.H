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
    Foam::CPCCellToCellExtStencil

Description

SourceFiles
    CPCCellToCellExtStencil.C

\*---------------------------------------------------------------------------*/

#ifndef CPCCellToCellExtStencil_H
#define CPCCellToCellExtStencil_H

#include "cellToCellExtStencil.H"
#include "boolList.H"
#include "HashSet.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class CPCCellToCellExtStencil Declaration
\*---------------------------------------------------------------------------*/

class CPCCellToCellExtStencil
:
    public cellToCellExtStencil
{
    // Private Member Functions

        //- Calculates per point the neighbour data (= pointCells)
        void calcPointBoundaryData
        (
            const boolList& isValidBFace,
            Map<labelList>& neiGlobal,
            Map<labelPairList>& neiTrafoGlobal
        ) const;

        void calcCellStencil(labelListList&, List<labelPairList>&) const;


public:

    // Constructors

        //- Construct from all cells and boundary faces
        explicit CPCCellToCellExtStencil(const polyMesh&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

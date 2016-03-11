/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
 2015-11-17 AtmosFOAM, Hilary Weller, University of Reading added support for
 cyclic boundaries
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
    Foam::cellToFaceExtStencil

Description
    Base class for extended cell-to-face stencils (face values from
    neighbouring cells)

SourceFiles
    cellToFaceExtStencil.C

\*---------------------------------------------------------------------------*/

#ifndef cellToFaceExtStencil_H
#define cellToFaceExtStencil_H

#include "globalIndex.H"
#include "globalIndexAndTransform.H"
#include "boolList.H"
#include "HashSet.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                      Class cellToFaceExtStencil Declaration
\*---------------------------------------------------------------------------*/

class cellToFaceExtStencil
{

public:

        //- Enumeration defining the location of the neighbour in the
        //  stencil
        enum neighbourLocation
        {
            NONE,               // N/A
            UNTRANSFORMED,      // after owner in untransformed part
            TRANSFORMED         // first element in transformed part
        };


private:

    // Private data

        const polyMesh& mesh_;

        //- Global numbering for cells and boundary faces
        const globalIndex globalNumbering_;


protected:

        //- Untransformed elements
        labelListList untransformedElements_;
       
        //- Transformed elements
        List<labelPairList> transformedElements_;


public:

        //- Merge two lists.
        //static void merge(const labelList&, labelList&);
        static void merge(const labelPairList&, labelPairList&);

        static void merge
        (
            const labelList& ownCCell,
            const labelList& neiCCells,

            labelList& faceStencil,
            labelHashSet& workSet
        );

        //- Helper: add transform to global cell labels
        static labelPairList transform
        (
            const polyMesh&,
            const label patchI,
            const globalIndex& globalNumbering,
            const labelList& untrafoStencil,
            const bool sendPatch
        );

        //- Helper: add transform to already transformed cells
        static labelPairList transform
        (
            const polyMesh&,
            const label patchI,
            const labelPairList& trafoStencil,
            const bool sendPatch
        );

        //- Collect cell neighbours into extended stencil
        void calcFaceStencil
        (
            const labelListList& untrafoCellCells,
            const List<labelPairList>& trafoCellCells,
            labelListList& faceStencil,
            List<labelPairList>& trafoFaceStencil
        ) const;

public:

    // Constructors

        //- Construct from mesh
        explicit cellToFaceExtStencil(const polyMesh&);


    // Member Functions

        const polyMesh& mesh() const
        {
            return mesh_;
        }

        //- Global numbering for cells and boundary faces
        const globalIndex& globalNumbering() const
        {
            return globalNumbering_;
        }

        //- Untransformed elements
        const labelListList& untransformedElements() const
        {
	        return untransformedElements_;
        }
        
        //- Untransformed elements
        labelListList& untransformedElements()
        {
            return untransformedElements_;
        }
        
        //- Transformed elements
        const List<labelPairList>& transformedElements() const
        {
            return transformedElements_;
        }
        
        //- Transformed elements
        List<labelPairList>& transformedElements()
        {
            return transformedElements_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

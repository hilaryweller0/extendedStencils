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
    Foam::cellToCellExtStencil

Description
    baseclass for extended cell centred addressing. Contains per cell a
    list of neighbouring cells and/or boundaryfaces in global addressing.

SourceFiles
    cellToCellExtStencil.C

\*---------------------------------------------------------------------------*/

#ifndef cellToCellExtStencil_H
#define cellToCellExtStencil_H

#include "globalIndex.H"
#include "boolList.H"
#include "HashSet.H"
#include "labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                      Class cellToCellExtStencil Declaration
\*---------------------------------------------------------------------------*/

class cellToCellExtStencil
{
    // Private data

        const polyMesh& mesh_;

        //- Global numbering for cells and boundary faces
        const globalIndex globalNumbering_;


protected:

        //- List (in globalNumbering) of elements without transformation
        labelListList untransformedElements_;

        //- List (in globalTransforms encoding) of elements with transformation
        List<labelPairList> transformedElements_;


        //- Merge two lists.
        static void merge
        (
            const label,
            const label,
            const labelList&,
            labelList&
        );

        //- Merge two lists.
        static void merge(const label, const labelList&, labelList&);

        //- Valid boundary faces (not empty and not coupled)
        void validBoundaryFaces(boolList& isValidBFace) const;

        //- Combine operator for labelLists
        template<class Type, class HashType>
        class unionEqOp
        {
            HashSet<Type, HashType>& work_;

            public:
                unionEqOp(HashSet<Type, HashType>& work)
                :
                    work_(work)
                {}
                void operator()(List<Type>& x, const List<Type>& y) const
                {
                    if (y.size())
                    {
                        if (x.empty())
                        {
                            x = y;
                        }
                        else
                        {
                            work_.clear();
                            work_.resize(2*(x.size()+y.size()));
                            work_.insert(x);
                            work_.insert(y);
                            x = work_.toc();
                        }
                    }
                };
        };

        //- Collect cell neighbours of faces in global numbering
        void insertFaceCells
        (
            const label exclude0,
            const label exclude1,
            const boolList& nonEmptyFace,
            const labelList& faceLabels,
            labelHashSet& globals
        ) const;

        //- Collect cell neighbours of faces in global numbering
        labelList calcFaceCells
        (
            const boolList& nonEmptyFace,
            const labelList& faceLabels,
            labelHashSet& globals
        ) const;


public:

    // Constructors

        //- Construct from mesh
        explicit cellToCellExtStencil(const polyMesh&);


    // Member Functions

        const labelListList& untransformedElements() const
        {
            return untransformedElements_;
        }

        const List<labelPairList>& transformedElements() const
        {
            return transformedElements_;
        }

        const polyMesh& mesh() const
        {
            return mesh_;
        }

        //- Global numbering for cells and boundary faces
        const globalIndex& globalNumbering() const
        {
            return globalNumbering_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

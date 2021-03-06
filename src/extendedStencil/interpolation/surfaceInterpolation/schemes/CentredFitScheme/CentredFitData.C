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

#include "CentredFitData.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "SVD.H"
#include "syncTools.H"
#include "extendedCentredCellToFaceExtStencil.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Polynomial>
Foam::CentredFitData<Polynomial>::CentredFitData
(
    const fvMesh& mesh,
    const extendedCentredCellToFaceExtStencil& stencil,
    const scalar linearLimitFactor,
    const scalar centralWeight
)
:
    FitData
    <
        CentredFitData<Polynomial>,
        extendedCentredCellToFaceExtStencil,
        Polynomial
    >
    (
        mesh, stencil, true, linearLimitFactor, centralWeight
    ),
    coeffs_(mesh.nFaces())
{
    if (debug)
    {
        InfoInFunction << "Contructing CentredFitData<Polynomial>" << endl;
    }

    calcFit();

    if (debug)
    {
        Info<<     "Finished constructing polynomialFit data" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Polynomial>
void Foam::CentredFitData<Polynomial>::calcFit()
{
    const fvMesh& mesh = this->mesh();

    const surfaceScalarField& w = mesh.surfaceInterpolation::weights();
    const surfaceScalarField::GeometricBoundaryField& bw = w.boundaryField();

    // Get the cell/face centres in stencil order.
    // Centred face stencils no good for triangles or tets.
    // Need bigger stencils
    List<List<point>> stencilPoints(mesh.nFaces());
    this->stencil().collectPositions(mesh.C(), stencilPoints);

    // find the fit coefficients for every face in the mesh

    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        FitData
        <
            CentredFitData<Polynomial>,
            extendedCentredCellToFaceExtStencil,
            Polynomial
        >::calcFit(coeffs_[facei], stencilPoints[facei], w[facei], facei);
    }

	// And for the boundaries
    forAll(bw, patchi)
    {
        const fvsPatchScalarField& pw = bw[patchi];

        if (pw.coupled())
        {
            label facei = pw.patch().start();

            forAll(pw, i)
            {
                FitData
                <
                    CentredFitData<Polynomial>,
                    extendedCentredCellToFaceExtStencil,
                    Polynomial
                >::calcFit(coeffs_[facei], stencilPoints[facei], pw[i], facei);
                facei++;
            }
        }
    }
}


// ************************************************************************* //

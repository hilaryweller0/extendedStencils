/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "UpwindFitData.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "SVD.H"
#include "syncTools.H"
#include "extendedUpwindCellToFaceExtStencil.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Polynomial>
Foam::UpwindFitData<Polynomial>::UpwindFitData
(
    const fvMesh& mesh,
    const extendedUpwindCellToFaceExtStencil& stencil,
    const bool linearCorrection,
    const scalar linearLimitFactor,
    const scalar centralWeight
)
:
    FitData
    <
        UpwindFitData<Polynomial>,
        extendedUpwindCellToFaceExtStencil,
        Polynomial
    >
    (
        mesh, stencil, linearCorrection, linearLimitFactor, centralWeight
    ),
    owncoeffs_(mesh.nFaces()),
    neicoeffs_(mesh.nFaces())
{
    if (debug)
    {
        Info<< "Contructing UpwindFitData<Polynomial>" << endl;
    }

    calcFit();

    if (debug)
    {
        Info<< "UpwindFitData<Polynomial>::UpwindFitData() :"
            << "Finished constructing polynomialFit data"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Polynomial>
void Foam::UpwindFitData<Polynomial>::calcFit()
{
    const fvMesh& mesh = this->mesh();

    const surfaceScalarField& w = mesh.surfaceInterpolation::weights();
    const surfaceScalarField::GeometricBoundaryField& bw = w.boundaryField();

    // Get the cell/face centres in stencil order.
    // Centred face stencils no good for triangles or tets.
    // Need bigger stencils
    List<List<point>> stencilPoints(mesh.nFaces());
    this->stencil().collectOwnPositions(mesh.C(), stencilPoints);

    // Owner stencil weights
    // ~~~~~~~~~~~~~~~~~~~~~

    // find the fit coefficients for every owner
    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        FitData
        <
            UpwindFitData<Polynomial>,
            extendedUpwindCellToFaceExtStencil,
            Polynomial
        >::calcFit(owncoeffs_[facei], stencilPoints[facei], w[facei], facei);
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
                    UpwindFitData<Polynomial>,
                    extendedUpwindCellToFaceExtStencil,
                    Polynomial
                >::calcFit
                (
                    owncoeffs_[facei], stencilPoints[facei], pw[i], facei
                );
                facei++;
            }
        }
    }


    // Neighbour stencil weights
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    // Note:reuse stencilPoints since is major storage
    this->stencil().collectNeiPositions(mesh.C(), stencilPoints);

    // find the fit coefficients for every neighbour
    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        FitData
        <
            UpwindFitData<Polynomial>,
            extendedUpwindCellToFaceExtStencil,
            Polynomial
        >::calcFit(neicoeffs_[facei], stencilPoints[facei], w[facei], facei);
    }

    forAll(bw, patchi)
    {
        const fvsPatchScalarField& pw = bw[patchi];

        if (pw.coupled())
        {
            label facei = pw.patch().start();

            forAll(pw, i)
            {
                scalarList wts(stencilPoints[facei].size(), scalar(1));
                FitData
                <
                    UpwindFitData<Polynomial>,
                    extendedUpwindCellToFaceExtStencil,
                    Polynomial
                >::calcFit
                (
                    neicoeffs_[facei], stencilPoints[facei], pw[i], facei
                );
                facei++;
            }
        }
    }
}


// ************************************************************************* //

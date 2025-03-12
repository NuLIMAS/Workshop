/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
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

#include "standingPressureWaveFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::standingPressureWaveFvPatchScalarField::
standingPressureWaveFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF)
{}


Foam::standingPressureWaveFvPatchScalarField::
standingPressureWaveFvPatchScalarField
(
    const standingPressureWaveFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
	A_(ptf.A_),
	k_(ptf.k_),
	T_(ptf.T_),
        off_(ptf.off_),
	lambda_(ptf.lambda_)
{}


Foam::standingPressureWaveFvPatchScalarField::
standingPressureWaveFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    A_(readScalar(dict.lookup("amplitude"))),
    k_(vector(dict.lookup("direction"))),
    T_(readScalar(dict.lookup("period"))),
    off_(readScalar(dict.lookup("offset"))),
    lambda_(readScalar(dict.lookup("wavelength")))
{}


Foam::standingPressureWaveFvPatchScalarField::
standingPressureWaveFvPatchScalarField
(
    const standingPressureWaveFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf),
	A_(ptf.A_),
	k_(ptf.k_),
	T_(ptf.T_),
    off_(ptf.off_),
	lambda_(ptf.lambda_)
{}


Foam::standingPressureWaveFvPatchScalarField::
standingPressureWaveFvPatchScalarField
(
    const standingPressureWaveFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF),
	A_(ptf.A_),
	k_(ptf.k_),
	T_(ptf.T_),
    off_(ptf.off_),
	lambda_(ptf.lambda_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::standingPressureWaveFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

	scalar pi = acos(-1.0);

	k_ = (2.0*pi/lambda_)*k_/mag(k_);

	const vectorField& x = patch().Cf();

	scalar omega = 2.0*pi/T_;

	const scalar t = db().time().value();

    operator==(off_+A_*cos(k_ & x) *cos(omega*t));

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void Foam::standingPressureWaveFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        standingPressureWaveFvPatchScalarField
    );
}

// ************************************************************************* //

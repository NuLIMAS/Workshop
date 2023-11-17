/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "propagatingPressureWaveFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::propagatingPressureWaveFvPatchScalarField::
propagatingPressureWaveFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF)
{}


Foam::propagatingPressureWaveFvPatchScalarField::
propagatingPressureWaveFvPatchScalarField
(
    const propagatingPressureWaveFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    k_(ptf.k_),
    T_(ptf.T_),
    wh_(ptf.wh_),
    wd_(ptf.wd_),    
    lambda_(ptf.lambda_),
    shift_(ptf.shift_),
    offset_(ptf.offset_)
{}


Foam::propagatingPressureWaveFvPatchScalarField::
propagatingPressureWaveFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    k_(vector(dict.lookup("direction"))),
    T_(readScalar(dict.lookup("period"))),
    wh_(readScalar(dict.lookup("waveHeight"))),
    wd_(readScalar(dict.lookup("waterDepth"))),
    lambda_(readScalar(dict.lookup("waveLength"))),
    shift_(dict.lookupOrDefault<scalar>("phase", 0.0)),
    offset_(dict.lookupOrDefault<scalar>("offset", 0.0))
{}


Foam::propagatingPressureWaveFvPatchScalarField::
propagatingPressureWaveFvPatchScalarField
(
    const propagatingPressureWaveFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf),
    k_(ptf.k_),
    T_(ptf.T_),
    wh_(ptf.wh_),
    wd_(ptf.wd_),
    lambda_(ptf.lambda_),
    shift_(ptf.shift_),
    offset_(ptf.offset_)
{}


Foam::propagatingPressureWaveFvPatchScalarField::
propagatingPressureWaveFvPatchScalarField
(
    const propagatingPressureWaveFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF),
    k_(ptf.k_),
    T_(ptf.T_),
    wh_(ptf.wh_),
    wd_(ptf.wd_),    
    lambda_(ptf.lambda_),
    shift_(ptf.shift_),
    offset_(ptf.offset_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::propagatingPressureWaveFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalar pi = acos(-1.0);

    const scalar t = db().time().value();

    k_ = (2.0*pi/lambda_)*k_/mag(k_);

    scalar A =0 ;

    A = 9810.0 * wh_ / (2.0* Foam::cosh(2.0*pi*wd_/lambda_));

    const vectorField& x = patch().Cf();

    scalar omega = 2.0*pi/T_;

    scalar th0 = 2.0*pi*shift_;

    operator==(offset_ + A*cos((k_ & x) - omega*t + th0));

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void Foam::propagatingPressureWaveFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("direction") << k_ << token::END_STATEMENT << nl;
    os.writeKeyword("waveHeight") << wh_ << token::END_STATEMENT << nl;
    os.writeKeyword("waterDepth") << wd_ << token::END_STATEMENT << nl;
    os.writeKeyword("waveLength") << lambda_ << token::END_STATEMENT << nl;
    os.writeKeyword("period") << T_ << token::END_STATEMENT << nl;
    os.writeKeyword("shift") <<shift_ << token::END_STATEMENT << nl;
    os.writeKeyword("offset") << offset_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        propagatingPressureWaveFvPatchScalarField
    );
}

// ************************************************************************* //

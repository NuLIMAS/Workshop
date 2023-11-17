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

#include "pressureWaveFixedValueFvPatchScalarField.H"
#include "solidContactFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pressureWaveFixedValueFvPatchScalarField::
pressureWaveFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF)
{}


Foam::pressureWaveFixedValueFvPatchScalarField::
pressureWaveFixedValueFvPatchScalarField
(
    const pressureWaveFixedValueFvPatchScalarField& ptf,
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


Foam::pressureWaveFixedValueFvPatchScalarField::
pressureWaveFixedValueFvPatchScalarField
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


Foam::pressureWaveFixedValueFvPatchScalarField::
pressureWaveFixedValueFvPatchScalarField
(
    const pressureWaveFixedValueFvPatchScalarField& ptf
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


Foam::pressureWaveFixedValueFvPatchScalarField::
pressureWaveFixedValueFvPatchScalarField
(
    const pressureWaveFixedValueFvPatchScalarField& ptf,
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

void Foam::pressureWaveFixedValueFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<vector>& U =
        patch().lookupPatchField<volVectorField, vector>("U");

    solidContactFvPatchVectorField& sc =
        const_cast<solidContactFvPatchVectorField&>
        (
            refCast<const solidContactFvPatchVectorField>(U)
        );
        


    vectorField masterPressure = sc.interpolateSlaveToMaster(sc.normalContactModelPtr()->slavePressure());
    
    scalar pi = acos(-1.0);

    k_ = (2.0*pi/lambda_)*k_/mag(k_);
    
    scalar A = 9810.0 * wh_ / (2.0* Foam::cosh(2.0*pi*wd_/lambda_));

    const vectorField& x = patch().Cf();

    scalar omega = 2.0*pi/T_;

    scalar th0 = 2.0*pi*shift_;

    const scalar t = db().time().value();

    //(offset_ + A*cos((k_ & x) - omega*t + th0))
    
  
    operator== ((offset_ + A*cos((k_ & x) - omega*t + th0)) + (patch().nf() & masterPressure)) ;

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void Foam::pressureWaveFixedValueFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    writeEntry("value", os);
    os.writeKeyword("waveHeight") << wh_ << token::END_STATEMENT << nl;
    os.writeKeyword("waterDepth") << wd_ << token::END_STATEMENT << nl;
    os.writeKeyword("direction") << k_ << token::END_STATEMENT << nl;
    os.writeKeyword("period") << T_ << token::END_STATEMENT << nl;
    os.writeKeyword("waveLength") << lambda_ << token::END_STATEMENT << nl;
    os.writeKeyword("phase") << shift_ << token::END_STATEMENT << nl;
    os.writeKeyword("offset") << offset_ << token::END_STATEMENT << nl;
    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        pressureWaveFixedValueFvPatchScalarField
    );
}

// ************************************************************************* //

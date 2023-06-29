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

#include "buoyantRhoPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

buoyantRhoPressureFvPatchScalarField::
buoyantRhoPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    rhoName_("rho")
{}


buoyantRhoPressureFvPatchScalarField::
buoyantRhoPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho"))
{
    fvPatchField<scalar>::operator=(patchInternalField());
    gradient() = 0.0;
}


buoyantRhoPressureFvPatchScalarField::
buoyantRhoPressureFvPatchScalarField
(
    const buoyantRhoPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    rhoName_(ptf.rhoName_)
{}


buoyantRhoPressureFvPatchScalarField::
buoyantRhoPressureFvPatchScalarField
(
    const buoyantRhoPressureFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    rhoName_(ptf.rhoName_)
{}


buoyantRhoPressureFvPatchScalarField::
buoyantRhoPressureFvPatchScalarField
(
    const buoyantRhoPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    rhoName_(ptf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void buoyantRhoPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // If variables are not found, evaluate as zero gradient
    // HJ, 17/Jan/2012
    /*if
    (
        //!db().foundObject<uniformDimensionedVectorField>("g")
        !db().foundObject<dimensionedVector>("g")
     || !db().foundObject<volScalarField>(rhoName_)
    )
    {
        InfoIn
        (
            "void buoyantRhoPressureFvPatchScalarField::updateCoeffs()"
        )   << "Fields required for evaluation not found for patch "
            << patch().name() << endl;

        gradient() = 0;
        fixedGradientFvPatchScalarField::updateCoeffs();

        return;
    }*/

    //const uniformDimensionedVectorField& g =
    //    db().lookupObject<uniformDimensionedVectorField>("g");
    
    const dictionary& materialProperties =
        db().lookupObject<IOdictionary>("materialProperties");
    dimensionedVector g(materialProperties.lookup("g"));

    const fvPatchField<scalar>& rho =
        lookupPatchField<volScalarField, scalar>("rho");
    const fvPatchField<scalar>& rhof =
        lookupPatchField<volScalarField, scalar>("rhof");
    const fvPatchField<scalar>& Cl =
        lookupPatchField<volScalarField, scalar>("Cl");
    const fvPatchField<scalar>& Cs =
        lookupPatchField<volScalarField, scalar>("Cs");
    //const fvPatchField<vector>& gradRho =
    //    lookupPatchField<volVectorField, vector>("gradRho");
    //const fvPatchField<scalar>& gh =
    //    lookupPatchField<volScalarField, scalar>("gh");

    // If the variable name is "p_rgh" or "pd" assume it is p - rho*g.h
    // and set the gradient appropriately.
    // Otherwise assume the variable is the static pressure.
    if
    (
        dimensionedInternalField().name() == "p_rgh"
     || dimensionedInternalField().name() == "pd"
    )
    {
        //gradient() = -rho*(g.value() & patch().nf());
        gradient() = (-rho.snGrad()*(g.value() & patch().Cf()) +  rhof.snGrad()*(g.value() & patch().Cf()))* Cl  ;
    }
    else
    {
        
        //gradient() = -rho*(g.value() & patch().nf());
        gradient() = (-rho*(g.value() & patch().Cf()) +  rhof*(g.value() & patch().Cf()))*Cl   ; //+ gh*(gradRho& patch().Cf())
    }

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void buoyantRhoPressureFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntryIfDifferent(os, "rho", word("rho"), rhoName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    buoyantRhoPressureFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

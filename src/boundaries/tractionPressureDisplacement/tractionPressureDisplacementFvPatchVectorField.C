/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
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

#include "tractionPressureDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tractionPressureDisplacementFvPatchVectorField::
tractionPressureDisplacementFvPatchVectorField
(
    const fvPatch& pr,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(pr, iF),
    traction_(pr.size(), vector::zero)
    //pressure_(p.size(), 0.0),
    //pName_("p")
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


tractionPressureDisplacementFvPatchVectorField::
tractionPressureDisplacementFvPatchVectorField
(
    const tractionPressureDisplacementFvPatchVectorField& tdpvf,
    const fvPatch& pr,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
   // const dictionary& dict
   
)
:
    fixedGradientFvPatchVectorField(tdpvf, pr, iF, mapper),
 //   pName_(dict.lookupOrDefault<word>("p", "p")),
    traction_(tdpvf.traction_, mapper)
    //pressure_(tdpvf.pressure_, mapper)
{}


tractionPressureDisplacementFvPatchVectorField::
tractionPressureDisplacementFvPatchVectorField
(
    const fvPatch& pr,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(pr, iF),
    traction_("traction", dict, pr.size())
   // pressure_("pressure", dict, p.size()),
   // pName_(tdpvf.pName_)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


tractionPressureDisplacementFvPatchVectorField::
tractionPressureDisplacementFvPatchVectorField
(
    const tractionPressureDisplacementFvPatchVectorField& tdpvf
)
:
    fixedGradientFvPatchVectorField(tdpvf),
    traction_(tdpvf.traction_)
    //pressure_(tdpvf.pressure_),
    //pName_(tdpvf.pName_)
{}


tractionPressureDisplacementFvPatchVectorField::
tractionPressureDisplacementFvPatchVectorField
(
    const tractionPressureDisplacementFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(tdpvf, iF),
    traction_(tdpvf.traction_)
    //pressure_(tdpvf.pressure_)
    //pName_(tdpvf.pName_)
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void tractionPressureDisplacementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
    traction_.autoMap(m);
    //pressure_.autoMap(m);
}


void tractionPressureDisplacementFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);

    const tractionPressureDisplacementFvPatchVectorField& dmptf =
        refCast<const tractionPressureDisplacementFvPatchVectorField>(ptf);

    traction_.rmap(dmptf.traction_, addr);
    //pressure_.rmap(dmptf.pressure_, addr);
}


void tractionPressureDisplacementFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    const dictionary& materialProperties =
        db().lookupObject<IOdictionary>("materialProperties");

    dimensionedScalar E(materialProperties.lookup("E"));
    dimensionedScalar nu(materialProperties.lookup("nu"));

    dimensionedScalar mu = E/(2.0*(1.0 + nu));
    dimensionedScalar lambda = nu*E/((1.0 + nu)*(1.0 - 2.0*nu));

    Switch planeStress(materialProperties.lookup("planeStress"));

    if (planeStress)
    {
        lambda = nu*E/((1.0 + nu)*(1.0 - nu));
    }

    scalar twoMuLambda = (2*mu + lambda).value();
       
    vectorField n = patch().nf(); //Unit normal

    const fvPatchField<symmTensor>& sigmaD =
        patch().lookupPatchField<volSymmTensorField, symmTensor>("sigmaD");

    const fvPatchField<scalar>& p =
        patch().lookupPatchField<volScalarField, scalar>("p");

    const volVectorField& U =
        patch().boundaryMesh().mesh().lookupObject<volVectorField>("U");

    vectorField pD =
        U.boundaryField()[patch().index()];

    gradient() =
    (
        (traction_ + (p)*n)
      + (twoMuLambda)*fvPatchField<vector>::snGrad() - (n & sigmaD)
    )/(twoMuLambda);




    fixedGradientFvPatchVectorField::updateCoeffs();
}


void tractionPressureDisplacementFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    traction_.writeEntry("traction", os);
   // pressure_.writeEntry("pressure", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    tractionPressureDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

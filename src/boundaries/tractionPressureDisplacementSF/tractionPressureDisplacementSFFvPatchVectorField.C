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

#include "tractionPressureDisplacementSFFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tractionPressureDisplacementSFFvPatchVectorField::
tractionPressureDisplacementSFFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


tractionPressureDisplacementSFFvPatchVectorField::
tractionPressureDisplacementSFFvPatchVectorField
(
    const tractionPressureDisplacementSFFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(tdpvf, p, iF, mapper),
    traction_(tdpvf.traction_, mapper),
    pressure_(tdpvf.pressure_, mapper)
{}


tractionPressureDisplacementSFFvPatchVectorField::
tractionPressureDisplacementSFFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_("traction", dict, p.size()),
    pressure_("pressure", dict, p.size())
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


tractionPressureDisplacementSFFvPatchVectorField::
tractionPressureDisplacementSFFvPatchVectorField
(
    const tractionPressureDisplacementSFFvPatchVectorField& tdpvf
)
:
    fixedGradientFvPatchVectorField(tdpvf),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_)
{}


tractionPressureDisplacementSFFvPatchVectorField::
tractionPressureDisplacementSFFvPatchVectorField
(
    const tractionPressureDisplacementSFFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(tdpvf, iF),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void tractionPressureDisplacementSFFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
    traction_.autoMap(m);
    pressure_.autoMap(m);
}


void tractionPressureDisplacementSFFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);

    const tractionPressureDisplacementSFFvPatchVectorField& dmptf =
        refCast<const tractionPressureDisplacementSFFvPatchVectorField>(ptf);

    traction_.rmap(dmptf.traction_, addr);
    pressure_.rmap(dmptf.pressure_, addr);
}


void tractionPressureDisplacementSFFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }



    vectorField n = patch().nf(); //Unit normal

    const fvPatchField<scalar>& mu =
        patch().lookupPatchField<volScalarField, scalar>("mu");
    const fvPatchField<scalar>& lambda =
        patch().lookupPatchField<volScalarField, scalar>("lambda");


    //const fvPatchField<symmTensor>& sigmaD =
    //    patch().lookupPatchField<volSymmTensorField, symmTensor>("sigmaD");

    //const fvPatchField<scalar>& p =
    //    patch().lookupPatchField<volScalarField, scalar>("p");

    const fvPatchField<tensor>& gradU =
        patch().lookupPatchField<volTensorField, tensor>("gradU");

    const dictionary& materialProperties =
        db().lookupObject<IOdictionary>("materialProperties");


    gradient() =
    (
        (traction_ + pressure_*n)
        - (n & (mu * gradU.T() - (mu + lambda)*gradU))
        - n * tr(gradU) * lambda
    )/ (2*mu + lambda);

    fixedGradientFvPatchVectorField::updateCoeffs();
}


void tractionPressureDisplacementSFFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeKeyword("patchType  splitCyclic;")<< nl;

    traction_.writeEntry("traction", os);
    pressure_.writeEntry("pressure", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    tractionPressureDisplacementSFFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

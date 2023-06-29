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


    const fvPatchField<symmTensor>& sigmaD =
        patch().lookupPatchField<volSymmTensorField, symmTensor>("sigmaD");

    const fvPatchField<scalar>& p =
        patch().lookupPatchField<volScalarField, scalar>("p");

    const fvPatchField<tensor>& gradU =
        patch().lookupPatchField<volTensorField, tensor>("gradU");
    //const fvPatchField<vector>& gradRho =
    //    patch().lookupPatchField<volVectorField, vector>("gradRho");
    //const fvPatchField<scalar>& gh =
    //    patch().lookupPatchField<volScalarField, scalar>("gh");
        
    //const fvPatchField<scalar>& rho =
    //    patch().lookupPatchField<volScalarField, scalar>("rho");
    //const fvPatchField<scalar>& rho =
    //    patch().lookupPatchField<volScalarField, scalar>("(((1-alpha)*rhof)+(alpha*rhos))");
    //const surfaceScalarField& rhodg =
    //    this->db().objectRegistry::lookupObject<surfaceScalarField>("rhodg");

    const dictionary& materialProperties =
        db().lookupObject<IOdictionary>("materialProperties");

    //const volScalarField& rhoref = 
    //    this->db().objectRegistry::lookupObject<volScalarField> ("rho");
    //const volVectorField& gradRhoref = 
    //    this->db().objectRegistry::lookupObject<volVectorField> ("gradRho");



    
    //scalarField  rholiq = rhoref.boundaryField()[5];
    //vectorField  gradRho = gradRhoref.boundaryField()[5];
    

    //dimensionedVector g(materialProperties.lookup("g"));
    //dimensionedScalar rhof(materialProperties.lookup("rhof"));
        

    //Info << "gh" << gh << endl;
    //Info << "gradRho" << gradRho << endl;
    //Info << "gradRho" << (gradRho*gh)   << endl;
    
    //Info << "patch cf" << g.component(vector::Z) * patch().Cf().component(vector::Z) << endl;
   
     
    //scalarField s1 =  patch().Cf().component(vector::Z)* gradRho.component(vector::Z)* -9.81;
    //Info << "s1" << s1 << endl;
    //vectorField  ref1 = (
    //    (traction_ + pressure_*n)
    //    - (n & (mu * gradU.T() - (mu + lambda)*gradU))
    //    - n * tr(gradU) * lambda
    //)/ (2*mu + lambda);
    
    //vectorField  ref2 = (
    //    (traction_ + pressure_*n)
    //  + (2*mu + lambda)*fvPatchField<vector>::snGrad() - (n & sigmaD)
    //  - p*n
    //)/(2*mu + lambda);
    //Info << "compare : " <<  ref1 << " " << ref2 << endl;
    //Info <<  mu << lambda << endl;
    
    //Info << mu<< endl;
    //Info << lambda<< endl;
    gradient() = 
    (
        (traction_ + pressure_*n)
        - (n & (mu * gradU.T() - (mu + lambda)*gradU))
        - n * tr(gradU) * lambda
        //+ (gh*gradRho) 
        //+ gradRho*gh
        //+ p*n
        //+ rho *(g.value() & patch().Cf())*n 
        //- n & rho*(g.value() & patch().nf())// pressure acting the opposite direction
    )/ (2*mu + lambda);
    //Info  << (rholiq-rho) * (g.value() & patch().Cf())*n  << endl;
    //Info  << (rho) << endl;
    
    //Info << "Patch " <<  patch().Cf() << endl;
    //Info << "rho" << rho << endl;
    //Info << "rhogh in tracPressure U" << rho*(g.value() & patch().Cf())*n << endl;
    
    //Info << "gradient" << gradient() << endl;
    //gradient() =
    //(
    //    (traction_ + pressure_*n)
    //  + (2*mu + lambda)*fvPatchField<vector>::snGrad() - (n & sigmaD)
    //)/(2*mu + lambda);

    // Pressure should be acting in opposite direction
    //gradient() -= (p*n)/(2*mu + lambda);

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

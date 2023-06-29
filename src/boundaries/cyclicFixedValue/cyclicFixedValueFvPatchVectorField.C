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
#include "splitCyclicFvPatchField.H"
#include "cyclicFixedValueFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"





// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::cyclicFixedValueFvPatchVectorField::
cyclicFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF)
{}


Foam::cyclicFixedValueFvPatchVectorField::
cyclicFixedValueFvPatchVectorField
(
    const cyclicFixedValueFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper)
{}


Foam::cyclicFixedValueFvPatchVectorField::
cyclicFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict)
{}


Foam::cyclicFixedValueFvPatchVectorField::
cyclicFixedValueFvPatchVectorField
(
    const cyclicFixedValueFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf)
{}


Foam::cyclicFixedValueFvPatchVectorField::
cyclicFixedValueFvPatchVectorField
(
    const cyclicFixedValueFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cyclicFixedValueFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    const fvPatchField<vector>& U =
        patch().lookupPatchField<volVectorField, vector>("U");
        
        
    //Info << "U vector : " << U << endl;
    
    //const fvMesh& mesh = patch().boundaryMesh().mesh();
    const volVectorField& Uref =this->db().objectRegistry::lookupObject<volVectorField>  ("U");
    
    
    vectorField  Uref2 = Uref.boundaryField()[6];

    vectorField n = patch().nf();

    // mesh.lookupObject<volVectorField>"U");

    //splitCyclicFvPatchField<vector>& sc =
    //    const_cast<splitCyclicFvPatchField<vector>&>
    //    (
    //        this->boundaryMesh()
    //    );
    
    //splitCyclicPolyPatch& sc =
    //    const_cast<splitCyclicPolyPatch&>
    //    (
    //        refCast<const splitCyclicPolyPatch>(U)
    //    );
    
    
    //Info << "label " << Ucyc << endl;
    
    //vectorField masterPressure = sc.interpolateSlaveToMaster(sc.normalContactModelPtr()->slavePressure());
         
     const vectorField x = patch().Cf(); 

    //vectorField Ucyc =
    //solidContactFvPatchVectorField& sc =
    //    const_cast<solidContactFvPatchVectorField&>
    //    (
    //        refCast<const solidContactFvPatchVectorField>(U)
    //    );
    //const fvPatchField<vector>& U =
    //    patch().lookupPatchField<volVectorField, vector>("U");

    //operator== Uref.boundaryField()[6];
    //operator == (Uref2);
    
    operator== (Uref2);

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::cyclicFixedValueFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchField<vector>::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        cyclicFixedValueFvPatchVectorField
    );
}

// ************************************************************************* //

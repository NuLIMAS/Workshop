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

    const volVectorField& Uref =this->db().objectRegistry::lookupObject<volVectorField>("U");


    int solBoundaryIndex = -1;
    for (size_t i = 0; i < Uref.boundaryField().size(); ++i)
    {
        if (Uref.boundaryField()[i].patch().name() == "sol")
        {
            solBoundaryIndex = i;
            break;
        }
    }
    if (solBoundaryIndex == -1)
    {
        Info << "Boundary field with name 'sol' not found." << endl;
        return;
    }
    vectorField  Uref2 = Uref.boundaryField()[solBoundaryIndex];

    vectorField n = patch().nf();

     const vectorField x = patch().Cf();


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

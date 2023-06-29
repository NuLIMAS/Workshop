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

#include "alphaInletOutletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "fvPatchField.H"
#include "fieldTypes.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::alphaInletOutletFvPatchScalarField::alphaInletOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    phiName_("phiAlpha")
{
    this->refValue() = pTraits<scalar>::zero;
    this->refGrad() = pTraits<scalar>::zero;
    this->valueFraction() = 0.0;
}



Foam::alphaInletOutletFvPatchScalarField::alphaInletOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF),
    phiName_(dict.lookupOrDefault<word>("phiAlpha", "phiAlpha"))
{
    // Read patch type
    this->readPatchType(dict);

    this->refValue() = Field<scalar>("inletValue", dict, p.size());



    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            Field<scalar>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(this->refValue());
    }

    this->refGrad() = pTraits<scalar>::zero;
    this->valueFraction() = 0.0;
}



Foam::alphaInletOutletFvPatchScalarField::alphaInletOutletFvPatchScalarField
(
    const alphaInletOutletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper),
    phiName_(ptf.phiName_)
{}



Foam::alphaInletOutletFvPatchScalarField::alphaInletOutletFvPatchScalarField
(
    const alphaInletOutletFvPatchScalarField& ptf
)
:
    mixedFvPatchField<scalar>(ptf),
    phiName_(ptf.phiName_)
{}



Foam::alphaInletOutletFvPatchScalarField::alphaInletOutletFvPatchScalarField
(
    const alphaInletOutletFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF),
    phiName_(ptf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::alphaInletOutletFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }


        const fvPatchField<scalar>& alpha =
            patch().lookupPatchField<volScalarField, scalar>("alpha");
        const scalarField& phiAlpha =
            patch().lookupPatchField<surfaceScalarField, scalar>("phiAlpha");
        const scalarField& alphaF =
            patch().lookupPatchField<surfaceScalarField, scalar>("alphaF");

        scalarField alpha_new = alpha;
        
        scalar s1 =0;
        scalar s2 =0;
        forAll(phiAlpha,facei)
        {
            if (phiAlpha[facei] < 0)
            {
                s1 += phiAlpha[facei];
                
            }
            if (phiAlpha[facei] > 0)
            {
                s2 += phiAlpha[facei]*alphaF[facei];
                //Info <<  "alpha " <<alphaF[facei] << endl;
            }

        }
        //Info << s2 << "  "<< s1 << endl;
        forAll(phiAlpha,facei)
        {
            if (phiAlpha[facei] < 0)
            {
                alpha_new[facei] = -s2/s1;
                this->refValue()[facei]=-s2/s1;
            }

        }


    // HR 2.1.19: Allow calls of virtual function in derived object registry
    // eg. postFixedSubRegistry
    if (!this->db().found(phiName_))
    {
        // Flux not available, do not update
        InfoIn
        (
            "void alphaInletOutletFvPatchField<Type>::"
            "updateCoeffs()"
        )   << "Flux field " << phiName_ << " not found.  "
            << "Performing mixed update" << endl;

        mixedFvPatchField<scalar>::updateCoeffs();

        return;
    }

    const scalarField& phip = this->lookupPatchField
    (
        phiName_,
        reinterpret_cast<const surfaceScalarField*>(0),
        reinterpret_cast<const scalar*>(0)
    );

    this->valueFraction() = 1.0 - pos(phip);

    mixedFvPatchField<scalar>::updateCoeffs();
}



void Foam::alphaInletOutletFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    if (phiName_ != "phiAlpha")
    {
        os.writeKeyword("phiAlpha")
            << phiName_ << token::END_STATEMENT << nl;
    }
    this->refValue().writeEntry("inletValue", os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

/*
void Foam::alphaInletOutletFvPatchScalarField::operator=
(
    const fvPatchField<scalar>& ptf
)
{
    fvPatchField<scalar>::operator=
    (
        this->valueFraction()*this->refValue()
        + (1 - this->valueFraction())*ptf
    );
}
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{    
    makePatchTypeField
    (
        fvPatchScalarField,
        alphaInletOutletFvPatchScalarField
    );

} // End namespace Foam

// ************************************************************************* //

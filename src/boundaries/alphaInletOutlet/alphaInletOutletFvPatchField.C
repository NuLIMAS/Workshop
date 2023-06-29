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

#include "alphaInletOutletFvPatchField.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
alphaInletOutletFvPatchField<Type>::alphaInletOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    phiName_("phiAlpha")
{
    this->refValue() = pTraits<Type>::zero;
    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 0.0;
}


template<class Type>
alphaInletOutletFvPatchField<Type>::alphaInletOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF),
    phiName_(dict.lookupOrDefault<word>("phiAlpha", "phiAlpha"))
{
    // Read patch type
    this->readPatchType(dict);

    this->refValue() = Field<Type>("inletValue", dict, p.size());



    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<Type>::operator=(this->refValue());
    }

    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 0.0;
}


template<class Type>
alphaInletOutletFvPatchField<Type>::alphaInletOutletFvPatchField
(
    const alphaInletOutletFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper),
    phiName_(ptf.phiName_)
{}


template<class Type>
alphaInletOutletFvPatchField<Type>::alphaInletOutletFvPatchField
(
    const alphaInletOutletFvPatchField<Type>& ptf
)
:
    mixedFvPatchField<Type>(ptf),
    phiName_(ptf.phiName_)
{}


template<class Type>
alphaInletOutletFvPatchField<Type>::alphaInletOutletFvPatchField
(
    const alphaInletOutletFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF),
    phiName_(ptf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void alphaInletOutletFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }


        //mixedFvPatchField<Type> alpha_new = this->refValue();
        
        // scalar s1 =0;
        // scalar s2 =0;
        // forAll(phiName_,facei)
        // {
        //     if (phiName_[facei] < 0)
        //     {
        //         //s1 += phiName_[facei];
                
        //     }
        //     if (phiName_[facei] > 0)
        //     {
        //         //s2 += phiName_[facei]*this->refValue()[facei];
        //         //Info <<  "alpha " <<alphaF[facei] << endl;
        //     }

        // }
        // //forAll(phiAlpha,facei)
        // //{
        // //    if (phiAlpha[facei] < 0)
        // //    {
        // //        alpha_new[facei] = -s2/s1;
        // //    }
        // //}


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

        mixedFvPatchField<Type>::updateCoeffs();

        return;
    }

    const scalarField& phip = this->lookupPatchField
    (
        phiName_,
        reinterpret_cast<const surfaceScalarField*>(0),
        reinterpret_cast<const scalar*>(0)
    );

    this->valueFraction() = 1.0 - pos(phip);

    mixedFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void alphaInletOutletFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    if (phiName_ != "phiAlpha")
    {
        os.writeKeyword("phiAlpha")
            << phiName_ << token::END_STATEMENT << nl;
    }
    this->refValue().writeEntry("inletValue", os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void alphaInletOutletFvPatchField<Type>::operator=
(
    const fvPatchField<Type>& ptf
)
{
    fvPatchField<Type>::operator=
    (
        this->valueFraction()*this->refValue()
        + (1 - this->valueFraction())*ptf
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

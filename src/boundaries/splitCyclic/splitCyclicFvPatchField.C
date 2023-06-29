/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
#include "transformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
splitCyclicFvPatchField<Type>::splitCyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    splitCyclicPatch_(refCast<const splitCyclicFvPatch>(p))
{}


template<class Type>
splitCyclicFvPatchField<Type>::splitCyclicFvPatchField
(
    const splitCyclicFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    splitCyclicPatch_(refCast<const splitCyclicFvPatch>(p))
{
    if (!isA<coupledFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "splitCyclicFvPatchField<Type>::splitCyclicFvPatchField"
            "("
                "const splitCyclicFvPatchField<Type>& ,"
                "const fvPatch&, "
                "const DimensionedField<Type, volMesh>&, "
                "const fvPatchFieldMapper&"
            ")"
        )   << "    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
splitCyclicFvPatchField<Type>::splitCyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<Type>(p, iF, dict),
    splitCyclicPatch_(refCast<const splitCyclicFvPatch>(p))
{
    if (!isA<splitCyclicFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "splitCyclicFvPatchField<Type>::splitCyclicFvPatchField"
            "("
                "const fvPatch&, "
                "const Field<Type>&, "
                "const dictionary&"
            ")",
            dict
        )   << "    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }

    this->evaluate(Pstream::blocking);
}


template<class Type>
splitCyclicFvPatchField<Type>::splitCyclicFvPatchField
(
    const splitCyclicFvPatchField<Type>& ptf
)
:
    splitCyclicLduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    splitCyclicPatch_(ptf.splitCyclicPatch_)
{}


template<class Type>
splitCyclicFvPatchField<Type>::splitCyclicFvPatchField
(
    const splitCyclicFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(ptf, iF),
    splitCyclicPatch_(ptf.splitCyclicPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type> > splitCyclicFvPatchField<Type>::patchNeighbourField() const
{
    const Field<Type>& iField = this->internalField();
    const labelUList& nbrFaceCells =
        splitCyclicPatch().splitCyclicPatch().neighbPatch().faceCells();

    tmp<Field<Type> > tpnf(new Field<Type>(this->size()));
    Field<Type>& pnf = tpnf();


    if (doTransform())
    {
        forAll(pnf, facei)
        {
            pnf[facei] = transform
            (
                forwardT()[0], iField[nbrFaceCells[facei]]
            );
        }
    }
    else
    {
        forAll(pnf, facei)
        {
            pnf[facei] = iField[nbrFaceCells[facei]];
        }
    }

    return tpnf;
}


template<class Type>
const splitCyclicFvPatchField<Type>& splitCyclicFvPatchField<Type>::neighbourPatchField()
const
{
    const GeometricField<Type, fvPatchField, volMesh>& fld =
    static_cast<const GeometricField<Type, fvPatchField, volMesh>&>
    (
        this->internalField()
    );

    return refCast<const splitCyclicFvPatchField<Type> >
    (
        fld.boundaryField()[this->coupledPatch().neighbPatchID()]
    );
}


template<class Type>
void splitCyclicFvPatchField<Type>::updateInterfaceMatrix
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes,
    const bool switchToLhs
) const
{
    const labelUList& nbrFaceCells =
        splitCyclicPatch().splitCyclicPatch().neighbPatch().faceCells();

    scalarField pnf(psiInternal, nbrFaceCells);

    // Transform according to the transformation tensors
    transformCoupleField(pnf, cmpt);

    // Multiply the field by coefficients and add into the result
    const labelUList& faceCells = splitCyclicPatch_.faceCells();

    forAll(faceCells, elemI)
    {
        result[faceCells[elemI]] -= coeffs[elemI]*pnf[elemI];
    }
}


template<class Type>
void Foam::splitCyclicFvPatchField<Type>::transformCoupleField
(
    scalarField& f,
    const direction cmpt
) const
{
    if (doTransform())
    {
        if (forwardT().size() == 1)
        {
            f *= pow(diag(forwardT()[0]).component(cmpt), rank());
        }
        else
        {
            f *= pow(diag(forwardT())().component(cmpt), rank());
        }
    }
}


template<class Type>
void splitCyclicFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

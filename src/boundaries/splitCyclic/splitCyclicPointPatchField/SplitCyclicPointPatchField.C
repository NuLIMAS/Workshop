/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "SplitCyclicPointPatchField.H"
#include "Swap.H"
#include "transformField.H"
#include "pointFields.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class SplitCyclicPointPatch,
    template<class> class MatrixType,
    class Type
>
SplitCyclicPointPatchField<PatchField, Mesh, PointPatch, SplitCyclicPointPatch, MatrixType, Type>::SplitCyclicPointPatchField
(
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF
)
:
    CoupledPointPatchField
    <
        PatchField,
        Mesh,
        PointPatch,
        typename SplitCyclicPointPatch::CoupledPointPatch,
        MatrixType,
        Type
    >(p, iF),
    splitCyclicPatch_(refCast<const SplitCyclicPointPatch>(p))
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class SplitCyclicPointPatch,
    template<class> class MatrixType,
    class Type
>
SplitCyclicPointPatchField<PatchField, Mesh, PointPatch, SplitCyclicPointPatch, MatrixType, Type>::SplitCyclicPointPatchField
(
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF,
    const dictionary& dict
)
:
    CoupledPointPatchField
    <
        PatchField,
        Mesh,
        PointPatch,
        typename SplitCyclicPointPatch::CoupledPointPatch,
        MatrixType,
        Type
    >(p, iF),
    splitCyclicPatch_(refCast<const SplitCyclicPointPatch>(p))
{
    if (!isType<SplitCyclicPointPatch>(p))
    {
        FatalIOErrorIn
        (
            "SplitCyclicPointPatchField<Type>::SplitCyclicPointPatchField\n"
            "(\n"
            "    const pointPatch& p,\n"
            "    const Field<Type>& field,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not SplitCyclic type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class SplitCyclicPointPatch,
    template<class> class MatrixType,
    class Type
>
SplitCyclicPointPatchField<PatchField, Mesh, PointPatch, SplitCyclicPointPatch, MatrixType, Type>::SplitCyclicPointPatchField
(
    const SplitCyclicPointPatchField
        <PatchField, Mesh, PointPatch, SplitCyclicPointPatch, MatrixType, Type>&,
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF,
    const PointPatchFieldMapper&
)
:
    CoupledPointPatchField
    <
        PatchField,
        Mesh,
        PointPatch,
        typename SplitCyclicPointPatch::CoupledPointPatch,
        MatrixType,
        Type
    >(p, iF),
    splitCyclicPatch_(refCast<const SplitCyclicPointPatch>(p))
{
    if (!isType<SplitCyclicPointPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "SplitCyclicPointPatchField<Type>::SplitCyclicPointPatchField\n"
            "(\n"
            "    const SplitCyclicPointPatchField<Type>& ptf,\n"
            "    const pointPatch& p,\n"
            "    const DimensionedField<Type, pointMesh>& iF,\n"
            "    const pointPatchFieldMapper& mapper\n"
            ")\n"
        )   << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class SplitCyclicPointPatch,
    template<class> class MatrixType,
    class Type
>
SplitCyclicPointPatchField<PatchField, Mesh, PointPatch, SplitCyclicPointPatch, MatrixType, Type>::SplitCyclicPointPatchField
(
    const SplitCyclicPointPatchField
    <PatchField, Mesh, PointPatch, SplitCyclicPointPatch, MatrixType, Type>& ptf
)
:
    CoupledPointPatchField
    <
        PatchField,
        Mesh,
        PointPatch,
        typename SplitCyclicPointPatch::CoupledPointPatch,
        MatrixType,
        Type
    >(ptf),
    splitCyclicPatch_(ptf.splitCyclicPatch_)
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class SplitCyclicPointPatch,
    template<class> class MatrixType,
    class Type
>
SplitCyclicPointPatchField
<PatchField, Mesh, PointPatch, SplitCyclicPointPatch, MatrixType, Type>::
SplitCyclicPointPatchField
(
    const SplitCyclicPointPatchField
    <PatchField, Mesh, PointPatch, SplitCyclicPointPatch, MatrixType, Type>& ptf,
    const DimensionedField<Type, Mesh>& iF
)
:
    CoupledPointPatchField
    <
        PatchField,
        Mesh,
        PointPatch,
        typename SplitCyclicPointPatch::CoupledPointPatch,
        MatrixType,
        Type
    >(ptf, iF),
    splitCyclicPatch_(ptf.splitCyclicPatch_)
{}

// ************************************************************************* //

}
// End namespace Foam

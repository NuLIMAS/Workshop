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

Class
    SplitCyclicPointPatchField

Description
    Cyclic front and back plane patch field

SourceFiles
    SplitCyclicPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef SplitCyclicPointPatchField_H
#define SplitCyclicPointPatchField_H

#include "CoupledPointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class SplitCyclicPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class SplitCyclicPointPatch,
    template<class> class MatrixType,
    class Type
>
class SplitCyclicPointPatchField
:
    public CoupledPointPatchField
    <
        PatchField,
        Mesh,
        PointPatch,
        typename SplitCyclicPointPatch::CoupledPointPatch,
        MatrixType,
        Type
    >
{
    // Private data

        //- Local reference cast into the SplitCyclic patch
        const SplitCyclicPointPatch& splitCyclicPatch_;


public:

    //- Runtime type information
    TypeName(SplitCyclicPointPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        SplitCyclicPointPatchField
        (
            const PointPatch&,
            const DimensionedField<Type, Mesh>&
        );

        //- Construct from patch, internal field and dictionary
        SplitCyclicPointPatchField
        (
            const PointPatch&,
            const DimensionedField<Type, Mesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        SplitCyclicPointPatchField
        (
            const SplitCyclicPointPatchField
            <PatchField, Mesh, PointPatch, SplitCyclicPointPatch, MatrixType, Type>&,
            const PointPatch&,
            const DimensionedField<Type, Mesh>&,
            const PointPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        SplitCyclicPointPatchField
        (
            const SplitCyclicPointPatchField
            <PatchField, Mesh, PointPatch, SplitCyclicPointPatch, MatrixType, Type>&
        );

        //- Construct and return a clone
        virtual autoPtr<PatchField<Type> > clone() const
        {
            return autoPtr<PatchField<Type> >
            (
                new SplitCyclicPointPatchField
                <
                    PatchField,
                    Mesh,
                    PointPatch,
                    SplitCyclicPointPatch,
                    MatrixType,
                    Type
                >(*this)
            );
        }

        //- Construct as copy setting internal field reference
        SplitCyclicPointPatchField
        (
            const SplitCyclicPointPatchField
            <PatchField, Mesh, PointPatch, SplitCyclicPointPatch, MatrixType, Type>&,
            const DimensionedField<Type, Mesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<PatchField<Type> > clone
        (
            const DimensionedField<Type, Mesh>& iF
        ) const
        {
            return autoPtr<PatchField<Type> >
            (
                new SplitCyclicPointPatchField
                <
                    PatchField,
                    Mesh,
                    PointPatch,
                    SplitCyclicPointPatch,
                    MatrixType,
                    Type
                >(*this, iF)
            );
        }


    // Member functions

        //- Constraint handling

            //- Return the constraint type this pointPatchField implements
            /*virtual const word& constraintType() const
            {
                return splitCyclicPointPatch::typeName;
            }*/

        //- Cyclic coupled interface functions

            //- Does the patch field perform the transfromation
            virtual bool doTransform() const
            {
                return !(splitCyclicPatch_.parallel() || pTraits<Type>::rank == 0);
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return splitCyclicPatch_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return splitCyclicPatch_.reverseT();
            }


        // Evaluation functions

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::blocking
            )
            {}

            //- Complete swap of patch point values and add to local values
            //virtual void swapAddSeparated
            //(
            //    const Pstream::commsTypes commsType,
            //    Field<Type>&
            //) const;

            //- Add field to internal field
            virtual void addField(Field<Type>&) const
            {}

        // Matrix construction (completion) functionality.  Currently faked
        // HJ, 31/May/2007

            //- Add diagonal/source contributions
            virtual void addDiag(scalarField&) const
            {}

            virtual void addSource(scalarField&) const
            {}

            //- Add upper/lower contributions
            virtual void addUpperLower(scalarField&) const
            {}

            virtual void eliminateUpperLower(scalarField&) const
            {}

            //- Get the cut edge coefficients in Amul order
            virtual tmp<scalarField> cutBouCoeffs(const lduMatrix&) const
            {
                return tmp<scalarField>(new scalarField(this->size(), 0));
            }

            //- Get the cut edge coefficients in Tmul order
            virtual tmp<scalarField> cutIntCoeffs(const lduMatrix&) const
            {
                return tmp<scalarField>(new scalarField(this->size(), 0));
            }

        // Coupled interface functionality

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const scalarField&,
                scalarField&,
                const lduMatrix&,
                const scalarField&,
                const direction,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const
            {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#    include "SplitCyclicPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

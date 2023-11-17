/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::processorSplitCyclicFvPatchField

Group
    grpCoupledBoundaryConditions

Description
    This boundary condition enables processor communication across cyclic
    patches.

    \heading Patch usage

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            processor;
    }
    \endverbatim

SeeAlso
    Foam::processorFvPatchField

SourceFiles
    processorSplitCyclicFvPatchField.C
    processorSplitCyclicFvPatchFields.H
    processorSplitCyclicFvPatchFields.C
    processorSplitCyclicFvPatchFieldsFwd.H

\*---------------------------------------------------------------------------*/

#ifndef processorSplitCyclicFvPatchField_H
#define processorSplitCyclicFvPatchField_H

#include "processorSplitCyclicFvPatch.H"
#include "processorFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class processorSplitCyclicFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class processorSplitCyclicFvPatchField
:
    public processorFvPatchField<Type>
{
    // Private data

        //- Local reference cast into the processor patch
        const processorSplitCyclicFvPatch& procPatch_;

public:

    //- Runtime type information
    TypeName(processorSplitCyclicFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        processorSplitCyclicFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch and internal field and patch field
        processorSplitCyclicFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        processorSplitCyclicFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given processorSplitCyclicFvPatchField onto a
        //  new patch
        processorSplitCyclicFvPatchField
        (
            const processorSplitCyclicFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        processorSplitCyclicFvPatchField(const processorSplitCyclicFvPatchField<Type>&);

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new processorSplitCyclicFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        processorSplitCyclicFvPatchField
        (
            const processorSplitCyclicFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new processorSplitCyclicFvPatchField<Type>(*this, iF)
            );
        }


    //- Destructor
    virtual ~processorSplitCyclicFvPatchField();


    // Member functions

        // Access

            //- Does the patch field perform the transfromation
            virtual bool doTransform() const
            {
                return !(procPatch_.parallel() || pTraits<Type>::rank == 0);
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return procPatch_.forwardT();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "processorSplitCyclicFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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

Class
    Foam::tractionDisplacementSFFvPatchVectorField

Description
    Fixed traction boundary condition for the standard linear elastic, fixed
    coefficient displacement equation.

SourceFiles
    tractionDisplacementSFFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef tractionDisplacementSFFvPatchVectorField_H
#define tractionDisplacementSFFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class tractionDisplacementFvPatch Declaration
\*---------------------------------------------------------------------------*/

class tractionDisplacementSFFvPatchVectorField
:
    public fixedGradientFvPatchVectorField
{

    // Private Data

        vectorField traction_;
        scalarField pressure_;


public:

    //- Runtime type information
    TypeName("tractionDisplacementSF");


    // Constructors

        //- Construct from patch and internal field
        tractionDisplacementSFFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        tractionDisplacementSFFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  tractionDisplacementSFFvPatchVectorField onto a new patch
        tractionDisplacementSFFvPatchVectorField
        (
            const tractionDisplacementSFFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        tractionDisplacementSFFvPatchVectorField
        (
            const tractionDisplacementSFFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new tractionDisplacementSFFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        tractionDisplacementSFFvPatchVectorField
        (
            const tractionDisplacementSFFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new tractionDisplacementSFFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            virtual const vectorField& traction() const
            {
                return traction_;
            }

            virtual vectorField& traction()
            {
                return traction_;
            }

            virtual const scalarField& pressure() const
            {
                return pressure_;
            }

            virtual  scalarField& pressure()
            {
                return pressure_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

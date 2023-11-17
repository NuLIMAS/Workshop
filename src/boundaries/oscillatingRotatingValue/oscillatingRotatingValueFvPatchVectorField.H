/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::oscillatingRotatingValueFvPatchField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides an oscillating condition in terms of
    amplitude and frequency.

        \f[
            x_p = (1 + a sin(2 \pi f t))x_{ref} + x_o
        \f]

    where

    \vartable
        x_p     | patch values
        x_{ref} | patch reference values
        x_o     | patch offset values
        a       | amplitude
        f       | frequency [1/s]
        t       | time [s]
    \endvartable

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        refValue     | reference value         | yes         |
        offset       | offset value            | no          | 0.0
        amplitude    | oscillation amplitude   | yes         |
        frequency    | oscillation frequency   | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            oscillatingRotatingValue;
        refValue        uniform 5.0;
        offset          0.0;
        amplitude       constant 0.5;
        frequency       constant 10;
    }
    \endverbatim

Note
    The amplitude and frequency entries are DataEntry types, able to describe
    time varying functions.  The example above gives the usage for supplying
    constant values.

SeeAlso
    Foam::DataEntry

SourceFiles
    oscillatingRotatingValueFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef oscillatingRotatingValueFvPatchVectorField_H
#define oscillatingRotatingValueFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class oscillatingRotatingValueFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/


class oscillatingRotatingValueFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Amplitude
        scalar amplitude_;

        //- Frequency
        scalar period_;
        
        //- rampPeriod
        scalar rampPeriod_;
        
        //- rampPeriod
        scalar width_;

        //- Centre of rotation
        vector origin_;

        //- Current time index
        label curTimeIndex_;

    // Private Member Functions

        //- Return current scale
        scalar currentScale() const;


public:

    //- Runtime type information
    TypeName("oscillatingRotatingValue");


    // Constructors

        //- Construct from patch and internal field
        oscillatingRotatingValueFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        oscillatingRotatingValueFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given oscillatingRotatingValueFvPatchField
        //  onto a new patch
        oscillatingRotatingValueFvPatchVectorField
        (
            const oscillatingRotatingValueFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        oscillatingRotatingValueFvPatchVectorField
        (
            const oscillatingRotatingValueFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new oscillatingRotatingValueFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        oscillatingRotatingValueFvPatchVectorField
        (
            const oscillatingRotatingValueFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new oscillatingRotatingValueFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return amplitude
            scalar amplitude() const
            {
                return amplitude_;
            }

            scalar& amplitude()
            {
                return amplitude_;
            }

            //- Return period
            scalar period() const
            {
                return period_;
            }

            scalar& period()
            {
                return period_;
            }

            //- Return rampPeriod
            scalar rampPeriod() const
            {
                return rampPeriod_;
            }

            scalar& rampPeriod()
            {
                return rampPeriod_;
            }
            
            //- Return plateWidth
            scalar width() const
            {
                return width_;
            }

            scalar& width()
            {
                return width_;
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
                const fvPatchField<vector>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

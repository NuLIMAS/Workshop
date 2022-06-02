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

Class
    Foam::propagatingPressureWaveFvPatchScalarField

Description
    Propagating pressure wave p = A*cos(k & x - omega*t + ).

    Example of the boundary condition specification:
    \verbatim
        inlet
        {
            type        propagatingPressureWave;
            amplitude    1.0;            // wave amplitude
            direction    (1 0 0);        // wave direction
            period        2.0;           // wave period
            wavelength    3.0;           // wave length
            phase        0.5;            // phase shift in rad/(2*pi)
            offset        0.0;           // base pressure
            value       uniform 0.0;     // initial value
        }
    \endverbatim

SourceFiles
    propagatingPressureWaveFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef propagatingPressureWaveFvPatchScalarField_H
#define propagatingPressureWaveFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class propagatingPressureWaveFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class propagatingPressureWaveFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Wave amplitude
        scalar A_;

        //- Wave direction
        vector k_;

        //- Wave period
        scalar T_;
        
        //- Wave height
        scalar wh_;
        
        //- Water depth
        scalar wd_;

        //- Wave length
        scalar lambda_;

        //- Wave phase shift
        scalar shift_;

        //- Mean pressure
        scalar offset_;


        

public:

    //- Runtime type information
    TypeName("propagatingPressureWave");


    // Constructors

        //- Construct from patch and internal field
        propagatingPressureWaveFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        propagatingPressureWaveFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  propagatingPressureWaveFvPatchScalarField
        //  onto a new patch
        propagatingPressureWaveFvPatchScalarField
        (
            const propagatingPressureWaveFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        propagatingPressureWaveFvPatchScalarField
        (
            const propagatingPressureWaveFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new propagatingPressureWaveFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        propagatingPressureWaveFvPatchScalarField
        (
            const propagatingPressureWaveFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new propagatingPressureWaveFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

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

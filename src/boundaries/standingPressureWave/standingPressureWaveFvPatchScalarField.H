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
    Foam::standingPressureWaveFvPatchScalarField

Description
    standing pressure wave p = A*cos(k&x)cos(omega*t)

    Example of the boundary condition specification:
    \verbatim
        inlet
        {
            type        standingPressureWave;
			amplitude	1.0;			// wave amplitude
			direction	(1 0 0);		// wave direction
			period		2.0;			// wave period
			wavelength	3.0;			// wave length
			phase		0.5;			// phase shift in rad/(2*pi)
			offset		0.0;			// base pressure
            value       uniform 0.0;    // initial value
        }
    \endverbatim

SourceFiles
    standingPressureWaveFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef standingPressureWaveFvPatchScalarField_H
#define standingPressureWaveFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class standingPressureWaveFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class standingPressureWaveFvPatchScalarField
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

        scalar off_;

		//- Wave length
		scalar lambda_;

public:

    //- Runtime type information
    TypeName("standingPressureWave");


    // Constructors

        //- Construct from patch and internal field
        standingPressureWaveFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        standingPressureWaveFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  standingPressureWaveFvPatchScalarField
        //  onto a new patch
        standingPressureWaveFvPatchScalarField
        (
            const standingPressureWaveFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        standingPressureWaveFvPatchScalarField
        (
            const standingPressureWaveFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new standingPressureWaveFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        standingPressureWaveFvPatchScalarField
        (
            const standingPressureWaveFvPatchScalarField&,
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
                new standingPressureWaveFvPatchScalarField
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

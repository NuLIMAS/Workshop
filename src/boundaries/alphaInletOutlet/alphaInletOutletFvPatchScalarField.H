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

Class
    Foam::alphaInletOutletFvPatchField

Description
    Foam::alphaInletOutletFvPatchField

SourceFiles
    alphaInletOutletFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef alphaInletOutletFvPatchScalarField_H
#define alphaInletOutletFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class alphaInletOutletFvPatch Declaration
\*---------------------------------------------------------------------------*/


class alphaInletOutletFvPatchScalarField
:
    public mixedFvPatchScalarField
{

protected:

    // Protected data

        word phiName_;


public:

    //- Runtime type information
    TypeName("alphaInletOutlet");


    // Constructors

        //- Construct from patch and internal field
        alphaInletOutletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        alphaInletOutletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given alphaInletOutletFvPatchField onto a new patch
        alphaInletOutletFvPatchScalarField
        (
            const alphaInletOutletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        alphaInletOutletFvPatchScalarField
        (
            const alphaInletOutletFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField > clone() const
        {
            return tmp<fvPatchScalarField >
            (
                new alphaInletOutletFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        alphaInletOutletFvPatchScalarField
        (
            const alphaInletOutletFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField > clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new alphaInletOutletFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the name of phi
            const word& phiName() const
            {
                return phiName_;
            }

            //- Return reference to the name of phi to allow adjustment
            word& phiName()
            {
                return phiName_;
            }


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

    //    virtual void operator=(const fvPatchScalarField& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

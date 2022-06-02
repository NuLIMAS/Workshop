/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "oscillatingRotatingValueFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "mathematicalConstants.H"
#include "septernion.H"
#include "quaternion.H"

using namespace Foam::mathematicalConstant;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


oscillatingRotatingValueFvPatchVectorField::oscillatingRotatingValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    amplitude_(),
    frequency_(),
    origin_(),
    curTimeIndex_(-1)
{}



oscillatingRotatingValueFvPatchVectorField::oscillatingRotatingValueFvPatchVectorField
(
    const oscillatingRotatingValueFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    amplitude_(ptf.amplitude_),
    frequency_(ptf.frequency_),  
    origin_(ptf.origin_),
    curTimeIndex_(-1)
{}



oscillatingRotatingValueFvPatchVectorField::oscillatingRotatingValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
    
)
:
    fixedValueFvPatchField<vector>(p, iF),
    amplitude_(readScalar(dict.lookup("amplitude"))),
    frequency_(readScalar(dict.lookup("frequency"))),
    origin_(vector(dict.lookup("origin"))),
    curTimeIndex_(-1)
{
    if (dict.found("value"))
    {
        fixedValueFvPatchField<vector>::operator==
        (
            Field<vector>("value", dict, p.size())
        );
    }
    else
    {
        const scalar t = db().time().value();

        // Convert the rotational motion from deg to rad
        scalar angle = asin(amplitude_/0.15);
        scalar eulerAngle = angle *sin(twoPi*frequency_* t);
        //eulerAngles *= pi/180.0;

        //quaternion R(eulerAngles.x(), eulerAngles.y(), eulerAngles.z());
        quaternion R(0 , eulerAngle, 0);
        septernion TR(septernion(origin_)*R*septernion(-origin_));
        
        vectorField displacement = this->patch().Cf();

        //Info << "displacement1  " << displacement << endl;
        forAll(displacement, faceI)
        {
            displacement[faceI] = TR.transform(displacement[faceI]);
        }
        //displacement = displacement - displacement1;
        displacement -= this->patch().Cf();
        //Info << "displacement 2 " << displacement << endl;

        fixedValueFvPatchField<vector>::operator== (displacement);
    }
}



oscillatingRotatingValueFvPatchVectorField::oscillatingRotatingValueFvPatchVectorField
(
    const oscillatingRotatingValueFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    amplitude_(ptf.amplitude_),
    frequency_(ptf.frequency_),  
    origin_(ptf.origin_),
    curTimeIndex_(-1)
{}



oscillatingRotatingValueFvPatchVectorField::oscillatingRotatingValueFvPatchVectorField
(
    const oscillatingRotatingValueFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    amplitude_(ptf.amplitude_),
    frequency_(ptf.frequency_),  
    origin_(ptf.origin_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void oscillatingRotatingValueFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<vector>::autoMap(m);
}



void oscillatingRotatingValueFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<vector>::rmap(ptf, addr);

    const oscillatingRotatingValueFvPatchVectorField& tiptf =
        refCast<const oscillatingRotatingValueFvPatchVectorField >(ptf);
}



void oscillatingRotatingValueFvPatchVectorField::updateCoeffs()


{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    
    {
    
        const scalar t = db().time().value();
        
        scalar angle = asin(amplitude_/0.15); //radian
        
        scalar eulerAngle = angle *sin(twoPi*frequency_* t);


        // Convert the rotational motion from deg to rad
        //eulerAngle *= pi/180.0;

        
        //quaternion R(eulerAngles.x(), eulerAngles.y(), eulerAngles.z());
        quaternion R(0,  eulerAngle, 0 );

        septernion TR(septernion(origin_)*R*septernion(-origin_));

        
        
        vectorField displacement = this->patch().Cf();

        forAll(displacement, faceI)
        {
            displacement[faceI] = TR.transform(displacement[faceI]);
        }

        displacement -= this->patch().Cf();
        fixedValueFvPatchField<vector>::operator== (displacement);


        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchField<vector>::updateCoeffs();
}



void oscillatingRotatingValueFvPatchVectorField::write(Ostream& os) const
{
    fixedValueFvPatchField<vector>::write(os);
    os.writeKeyword("amplitude") << amplitude_ << token::END_STATEMENT << nl;
    os.writeKeyword("frequency") << frequency_ << token::END_STATEMENT << nl;
    os.writeKeyword("origin") << origin_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        oscillatingRotatingValueFvPatchVectorField
    );
}

// ************************************************************************* //

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
    Foam::splitCyclicLduInterface

Description
    An abstract base class for splitCyclic coupled interfaces.

SourceFiles
    splitCyclicLduInterface.C

\*---------------------------------------------------------------------------*/

#ifndef splitCyclicLduInterface_H
#define splitCyclicLduInterface_H

#include "lduInterface.H"
#include "primitiveFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class splitCyclicFvPatch;

/*---------------------------------------------------------------------------*\
                  Class splitCyclicLduInterface Declaration
\*---------------------------------------------------------------------------*/

class splitCyclicLduInterface
{

public:

    //- Runtime type information
    TypeName("splitCyclicLduInterface");

    // Constructors

    //- Destructor
    virtual ~splitCyclicLduInterface();


    // Member Functions

        // Access

            //- Neighbour patchID
            virtual label neighbPatchID() const = 0;

            //- Neighbour patchID
            virtual bool master() const = 0;

            //virtual const splitCyclicLduInterface& neighbInterface() const = 0;

            virtual const splitCyclicLduInterface& neighbPatch() const = 0;

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const = 0;

            //- Return face reverse transformation tensor
            virtual const tensorField& reverseT() const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

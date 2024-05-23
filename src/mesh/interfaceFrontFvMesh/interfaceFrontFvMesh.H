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
    interfaceFrontFvMesh

Description
    Splits the internal face between liqufied and non-liqufied regions

SourceFiles
    interfaceFrontFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceFrontFvMesh_H
#define interfaceFrontFvMesh_H


#include "topoChangerFvMesh.H"
#include "interfaceFront.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class interfaceFrontFvMesh Declaration
\*---------------------------------------------------------------------------*/


class interfaceFrontFvMesh
:
    public topoChangerFvMesh
{
    // Private data

        //- Motion dictionary
        dictionary motionDict_;

        //- Face zone ID
        label openFzID_;

        //- Face zone ID
        label wallFzID_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        interfaceFrontFvMesh(const interfaceFrontFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const interfaceFrontFvMesh&);


        //- Add attach/detach zones and modifiers
        void addZonesAndModifiers();

        //- correct boundary patches
        void correctPatches();


public:
     //- Runtime type information
         TypeName("interfaceFrontFvMesh");


     // Constructors

        //- Construct from database
        explicit interfaceFrontFvMesh(const IOobject& io);


    // Destructor

        virtual ~interfaceFrontFvMesh();


    // Access Functions

        const dictionary& motionDict() const
        {
            return motionDict_;
        }


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

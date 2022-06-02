/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    simpleLiqFvMesh

Description

SourceFiles
    simpleLiqFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef simpleLiqFvMesh_H
#define simpleLiqFvMesh_H

 
#include "topoChangerFvMesh.H"
#include "attachDetach.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class simpleLiqFvMesh Declaration
\*---------------------------------------------------------------------------*/


class simpleLiqFvMesh
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
        simpleLiqFvMesh(const simpleLiqFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const simpleLiqFvMesh&);


        //- Add attach/detach zones and modifiers
        void addZonesAndModifiers();

        //- correct boundary patches
        void correctPatches();


public:
        //- Runtime type information
         TypeName("simpleLiqFvMesh");
   
   
     // Constructors

        //- Construct from database
        explicit simpleLiqFvMesh(const IOobject& io);


    // Destructor

        virtual ~simpleLiqFvMesh();


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

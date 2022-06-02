/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::cyclicPointPatch

Description
    Cyclic patch for post-processing.

SourceFiles
    splitCyclicPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef splitCyclicPointPatch_H
#define splitCyclicPointPatch_H

#include "coupledFacePointPatch.H"
#include "splitCyclicPolyPatch.H"
#include "pointBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class cyclicPointPatch Declaration
\*---------------------------------------------------------------------------*/

class splitCyclicPointPatch
:
    public coupledFacePointPatch
{
    // Private data

        //- Local reference cast into the cyclic patch
        const splitCyclicPolyPatch& splitCyclicPolyPatch_;


    // Private Member Functions

        //- Disallow default construct as copy
        splitCyclicPointPatch(const splitCyclicPointPatch&);

        //- Disallow default assignment
        void operator=(const splitCyclicPointPatch&);


    // Demand driven private data

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry();

        //- Calculate the patch geometry
        virtual void calcGeometry();

        //- Initialise the patches for moving pointsq
        virtual void initMovePoints(const pointField&);

        //- Correct patches after moving points
        virtual void movePoints( );

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh();

        //- Update of the patch topology
        virtual void updateMesh();


public:

    //- Runtime type information
    TypeName(splitCyclicPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        splitCyclicPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        );


    //- Destructor
    virtual ~splitCyclicPointPatch();


    // Member Functions

        // Access

            //- Return the constraint type this pointPatch implements.
            virtual const word& constraintType() const
            {
                return type();
            }

            //- Return the underlying splitCyclicPolyPatch
            const splitCyclicPolyPatch& splitCyclicPatch() const
            {
                return splitCyclicPolyPatch_;
            }

            //- Return neighbour point patch
            const splitCyclicPointPatch& neighbPatch() const
            {
                label patchI = splitCyclicPolyPatch_.neighbPatchID();
                const pointPatch& pp = this->boundaryMesh()[patchI];
                return refCast<const splitCyclicPointPatch>(pp);
            }

            //- Are the cyclic planes parallel
            bool parallel() const
            {
                return splitCyclicPolyPatch_.parallel();
            }

            //- Return face transformation tensor
            const tensorField& forwardT() const
            {
                return splitCyclicPolyPatch_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            const tensorField& reverseT() const
            {
                return splitCyclicPolyPatch_.reverseT();
            }


        // Access functions for demand driven data

            //- Return the set of pairs of points that require transformation
            //  and/or mapping. First index is on this patch, second on the
            //  neighbour patch.
            virtual const edgeList& transformPairs() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::processorSplitCyclicFvPatch

Description
    Processor patch.

SourceFiles
    processorSplitCyclicFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef processorSplitCyclicFvPatch_H
#define processorSplitCyclicFvPatch_H

#include "processorSplitCyclicPolyPatch.H"
#include "processorFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class processorSplitCyclicFvPatch Declaration
\*---------------------------------------------------------------------------*/

class processorSplitCyclicFvPatch
:
    public processorFvPatch
{
    // Private Data

        const processorSplitCyclicPolyPatch& procPolyPatch_;

public:

    //- Runtime type information
    TypeName(processorSplitCyclicPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        processorSplitCyclicFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            processorFvPatch(patch, bm),
            procPolyPatch_(refCast<const processorSplitCyclicPolyPatch>(patch))
        {}


    // Member functions

        //- Return message tag used for sending
        virtual int tag() const
        {
            // Allocate from Pstream?
            return procPolyPatch_.tag();
        }

        const processorSplitCyclicPolyPatch& procPolyPatch() const
        {
            return procPolyPatch_;
        }

        //- Are the cyclic planes parallel
        virtual bool parallel() const
        {
            return procPolyPatch_.parallel();
        }

        //- Return face transformation tensor
        virtual const tensorField& forwardT() const
        {
            return procPolyPatch_.forwardT();
        }

        //- Return neighbour-cell transformation tensor
        virtual const tensorField& reverseT() const
        {
            return procPolyPatch_.reverseT();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

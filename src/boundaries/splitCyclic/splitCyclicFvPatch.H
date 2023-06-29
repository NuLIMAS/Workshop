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
    Foam::splitCyclicFvPatch

Description
    splitCyclic-plane patch.

SourceFiles
    splitCyclicFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef splitCyclicFvPatch_H
#define splitCyclicFvPatch_H

#include "coupledFvPatch.H"
#include "splitCyclicLduInterface.H"
#include "splitCyclicPolyPatch.H"
#include "fvBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class splitCyclicFvPatch Declaration
\*---------------------------------------------------------------------------*/

class splitCyclicFvPatch
:
    public coupledFvPatch,
    public splitCyclicLduInterface
{
    // Private data

        const splitCyclicPolyPatch& splitCyclicPolyPatch_;


protected:

    // Protected Member functions

        //- Make patch weighting factors
        void makeWeights(fvsPatchScalarField&) const;

        //- Make patch face - neighbour cell distances
        void makeDeltaCoeffs(fvsPatchScalarField&) const;

public:

    //- Runtime type information
    TypeName(splitCyclicPolyPatch::typeName_());


    // Constructors

        //- Construct from splitPolyPatch
        splitCyclicFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            coupledFvPatch(patch, bm),
            splitCyclicPolyPatch_(refCast<const splitCyclicPolyPatch>(patch))
        {}


    // Member functions

        // Access

            //- Return local reference cast into the splitCyclic patch
            const splitCyclicPolyPatch& splitCyclicPatch() const
            {
                return splitCyclicPolyPatch_;
            }

            //- Return neighbour
            virtual label neighbPatchID() const
            {
                return splitCyclicPolyPatch_.neighbPatchID();
            }

            virtual bool master() const
            {
                return splitCyclicPolyPatch_.master();
            }

            //- Return processor number
            virtual const splitCyclicFvPatch& neighbPatch() const
            {
                Info << this->boundaryMesh()[splitCyclicPolyPatch_.neighbPatchID()].type() << endl;
                Info << typeid(this->boundaryMesh()[splitCyclicPolyPatch_.neighbPatchID()]).name() << endl;
                return refCast<const splitCyclicFvPatch>
                (
                    this->boundaryMesh()[splitCyclicPolyPatch_.neighbPatchID()]
                );
                //return reinterpret_cast<const splitCyclicFvPatch&>
                //(
                //    this->boundaryMesh()[splitCyclicPolyPatch_.neighbPatchID()]
                //);
                //return this->boundaryMesh()[splitCyclicPolyPatch_.neighbPatchID()];
            }

            //- Are the splitCyclic planes parallel
            virtual bool parallel() const
            {
                return splitCyclicPolyPatch_.parallel();
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return splitCyclicPolyPatch_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return splitCyclicPolyPatch_.reverseT();
            }

            virtual const splitCyclicFvPatch& neighbFvPatch() const
            {
                return refCast<const splitCyclicFvPatch>
                (
                    this->boundaryMesh()[splitCyclicPolyPatch_.neighbPatchID()]
                );
            }


            //- Return delta (P to N) vectors across coupled patch
            virtual tmp<vectorField> delta() const;


        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const;

            //- Transfer and return neighbour field
            virtual tmp<labelField> transfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& interfaceData
            ) const;

            //- Return neighbour field
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& internalData
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

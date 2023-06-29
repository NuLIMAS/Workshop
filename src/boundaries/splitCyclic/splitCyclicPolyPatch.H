/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::splitCyclicPolyPatch

Description
    splitCyclic plane patch.

    Note: morph patch face ordering uses geometric matching so with the
    following restrictions:
        -coupled patches should be flat planes.
        -no rotation in patch plane

    Uses splitCoupledPolyPatch::calcFaceTol to calculate
    tolerance per face which might need tweaking.

    Switch on 'splitCyclicPolyPatch' debug flag to write .obj files to show
    the matching.

SourceFiles
    splitCyclicPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef splitCyclicPolyPatch_H
#define splitCyclicPolyPatch_H

#include "coupledPolyPatch.H"
#include "edgeList.H"
#include "polyBoundaryMesh.H"
#include "diagTensorField.H"
#include "coupleGroupIdentifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class splitCyclicPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class splitCyclicPolyPatch
:
    public coupledPolyPatch
{
public:

    enum transformType
    {
        UNKNOWN,
        ROTATIONAL,
        TRANSLATIONAL,
        COINCIDENTFULLMATCH,// assume no transforms
                            // and check the points in faces match
        NOORDERING          // unspecified, no automatic ordering

    };
    static const NamedEnum<transformType, 5> transformTypeNames;

private:

    // Private data

        //- Name of other half
        mutable word neighbPatchName_;

        //- Optional patchGroup to find neighbPatch
        const coupleGroupIdentifier coupleGroup_;

        //- Index of other half
        mutable label neighbPatchID_;

        //- Type of transformation - rotational or translational
        transformType transform_;

        // For rotation

            //- Axis of rotation for rotational splitCyclics
            vector rotationAxis_;

            //- Point on axis of rotation for rotational splitCyclics
            point rotationCentre_;

        // For translation

            //- Translation vector
            vector separationVector_;


        //- List of edges formed from connected points. e[0] is the point on
        //  the first half of the patch, e[1] the corresponding point on the
        //  second half.
        mutable edgeList* coupledPointsPtr_;

        //- List of connected edges. e[0] is the edge on the first half of the
        //  patch, e[1] the corresponding edge on the second half.
        mutable edgeList* coupledEdgesPtr_;

        //- Temporary storage of owner side patch during ordering.
        mutable autoPtr<primitivePatch> ownerPatchPtr_;


    // Private Member Functions

        //- Find amongst selected faces the one with the largest area
        static label findMaxArea(const pointField&, const faceList&);

        void calcTransforms
        (
            const primitivePatch& half0,
            const pointField& half0Ctrs,
            const vectorField& half0Areas,
            const pointField& half1Ctrs,
            const vectorField& half1Areas
        );

        // Face ordering

            // Given a split of faces into left and right half calculate the
            // centres and anchor points. Transform the left points so they
            // align with the right ones
            void getCentresAndAnchors
            (
                const primitivePatch& pp0,
                const primitivePatch& pp1,

                pointField& half0Ctrs,
                pointField& half1Ctrs,
                pointField& anchors0,
                scalarField& tols
            ) const;

            //- Return normal of face at max distance from rotation axis
            vector findFaceMaxRadius(const pointField& faceCentres) const;


protected:

    // Protected Member functions

        //- Initialise the calculation of the patch addressing
        virtual void initAddressing();

        //- Calculate the patch addressing
        virtual void calcAddressing();

        //- Recalculate the transformation tensors
        virtual void calcTransforms();

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry();

        //- Initialise the calculation of the patch geometry
        //virtual void initGeometry
        //(
        //    const primitivePatch& referPatch,
        //    pointField& nbrCtrs,
        //    vectorField& nbrAreas,
        //    pointField& nbrCc
        //);

        //- Calculate the patch geometry
        virtual void calcGeometry();

        //- Calculate the patch geometry
        // virtual void calcGeometry
        // (
        //     const primitivePatch& referPatch,
        //     const pointField& thisCtrs,
        //     const vectorField& thisAreas,
        //     const pointField& thisCc,
        //     const pointField& nbrCtrs,
        //     const vectorField& nbrAreas,
        //     const pointField& nbrCc
        // );

        //- Initialise the patches for moving points
        virtual void initMovePoints(const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh();

        //- Update of the patch topology
        virtual void updateMesh();

public:

    //- Declare friendship with processorsplitCyclicPolyPatch
    friend class processorsplitCyclicPolyPatch;


    //- Runtime type information
    TypeName("splitCyclic");


    // Constructors

        //- Construct from components
        splitCyclicPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm
        );

        //- Construct from dictionary
        splitCyclicPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm
        );

        //- Construct as copy
        splitCyclicPolyPatch(const splitCyclicPolyPatch&);

        //- Construct as copy, resetting the boundary mesh
        splitCyclicPolyPatch(const splitCyclicPolyPatch&, const polyBoundaryMesh&);

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        splitCyclicPolyPatch
        (
            const splitCyclicPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );


    //- Destructor
    virtual ~splitCyclicPolyPatch();


    // Member Functions

        //- Neighbour patch name
        const word& neighbPatchName() const;

        //- Neighbour patchID
        label neighbPatchID() const;

        //- Is this the master side?
        //virtual bool master() const
        //{
        //    notImplemented("missing functionality");
        //    return true;
        //}
        virtual bool master() const
        {
            return index() < neighbPatchID();
        }

        virtual bool slave() const
        {
            return !master();
        }

        const splitCyclicPolyPatch& neighbPatch() const
        {
            const polyPatch& pp = this->boundaryMesh()[neighbPatchID()];
            return refCast<const splitCyclicPolyPatch>(pp);
        }

        //- Return connected points (from patch local to neighbour patch local)
        //  Demand driven calculation. Does primitivePatch::clearOut after
        //  calculation!
        const edgeList& coupledPoints() const;

        //- Return connected edges (from patch local to neighbour patch local).
        //  Demand driven calculation. Does primitivePatch::clearOut after
        //  calculation!
        const edgeList& coupledEdges() const;

        //- Transform a patch-based position from other side to this side
        //virtual void transformPosition(pointField& l) const;

        //- Transform a patch-based position from other side to this side
        //virtual void transformPosition(point&, const label facei) const;


        // Transformation

        label transformGlobalFace(const label facei) const
        {
            label offset = facei-start();
            label neighbStart = neighbPatch().start();

            if (offset >= 0 && offset < size())
            {
                return neighbStart+offset;
            }
            else
            {
                FatalErrorIn
                (
                    "splitCyclicPolyPatch::transformGlobalFace(const label) const"
                )   << "Face " << facei << " not in patch " << name()
                    << exit(FatalError);
                return -1;
            }
        }

        //- Type of transform
        transformType transform() const
        {
            return transform_;
        }

        //- Axis of rotation for rotational splitCyclics
        const vector& rotationAxis() const
        {
            return rotationAxis_;
        }

        //- Point on axis of rotation for rotational splitCyclics
        const point& rotationCentre() const
        {
            return rotationCentre_;
        }

        //- Translation vector for translational splitCyclics
        const vector& separationVector() const
        {
            return separationVector_;
        }


        //- Initialize ordering for primitivePatch. Does not
        //  refer to *this (except for name() and type() etc.)
        virtual void initOrder( const primitivePatch&) const;

        //- Return new ordering for primitivePatch.
        //  Ordering is -faceMap: for every face
        //  index of the new face -rotation:for every new face the clockwise
        //  shift of the original face. Return false if nothing changes
        //  (faceMap is identity, rotation is 0), true otherwise.
        virtual bool order
        (
            const primitivePatch&,
            labelList& faceMap,
            labelList& rotation
        ) const;

        //- Synchronise communications of ordering for primitivePatch
        //  Used in cases when no topological change happens locally,
        //  but is happening on other processors
        virtual void syncOrder() const;

        //- Write the splitPolyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

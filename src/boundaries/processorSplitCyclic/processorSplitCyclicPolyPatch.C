/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "processorSplitCyclicPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "SubField.H"
#include "splitCyclicPolyPatch.H"
#include "dictionary.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(processorSplitCyclicPolyPatch, 0);
    addToRunTimeSelectionTable(polyPatch, processorSplitCyclicPolyPatch, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorSplitCyclicPolyPatch::processorSplitCyclicPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const int myProcNo,
    const int neighbProcNo
    //const word& referPatchName,
    //const transformType transform,
    //const word& patchType
)
:
    processorPolyPatch
    (
        name,
        size,
        start,
        index,
        bm,
        myProcNo,
        neighbProcNo
        //transform,
        //patchType
    ),
    referPatchName_(referPatchName()),
    tag_(-1),
    referPatchID_(-1)
{}


Foam::processorSplitCyclicPolyPatch::processorSplitCyclicPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
    //const word& patchType
)
:
    processorPolyPatch(name, dict, index, bm),
    referPatchName_(dict.lookup("referPatch")),
    tag_(dict.lookupOrDefault<int>("tag", -1)),
    referPatchID_(-1)
{}


Foam::processorSplitCyclicPolyPatch::processorSplitCyclicPolyPatch
(
    const processorSplitCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    processorPolyPatch(pp, bm),
    referPatchName_(pp.referPatchName()),
    tag_(pp.tag()),
    referPatchID_(-1)
{}


Foam::processorSplitCyclicPolyPatch::processorSplitCyclicPolyPatch
(
    const processorSplitCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    processorPolyPatch(pp, bm, index, newSize, newStart),
    referPatchName_(pp.referPatchName_),
    tag_(pp.tag()),
    referPatchID_(-1)
{}


Foam::processorSplitCyclicPolyPatch::processorSplitCyclicPolyPatch
(
    const processorSplitCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& referPatchName
)
:
    processorPolyPatch(pp, bm, index, newSize, newStart),
    referPatchName_(referPatchName),
    tag_(-1),
    referPatchID_(-1)
{}


/*Foam::processorSplitCyclicPolyPatch::processorSplitCyclicPolyPatch
(
    const processorSplitCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    processorPolyPatch(pp, bm, index, mapAddressing, newStart),
    referPatchName_(pp.referPatchName()),
    tag_(-1),
    referPatchID_(-1)
{}*/


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::processorSplitCyclicPolyPatch::~processorSplitCyclicPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

int Foam::processorSplitCyclicPolyPatch::tag() const
{
    if (tag_ == -1)
    {
        // Get unique tag to use for all comms. Make sure that both sides
        // use the same tag
        const splitCyclicPolyPatch& cycPatch = refCast<const splitCyclicPolyPatch>
        (
            referPatch()
        );
        if (master())
        {
            tag_ = Hash<word>()(cycPatch.name()) % 32768u;
        }
        else
        {
            tag_ = Hash<word>()(cycPatch.neighbPatch().name()) % 32768u;
        }

        if (tag_ == Pstream::msgType() || tag_ == -1)
        {
            FatalErrorIn("processorSplitCyclicPolyPatch::tag() const")
                << "Tag calculated from cyclic patch name " << tag_
                << " is the same as the current message type "
                << Pstream::msgType() << " or -1" << nl
                << "Please set a non-conflicting, unique, tag by hand"
                << " using the 'tag' entry"
                << exit(FatalError);
        }
        if (debug)
        {
            Pout<< "processorSplitCyclicPolyPatch " << name() << " uses tag " << tag_
                << endl;
        }
    }
    return tag_;
}


void Foam::processorSplitCyclicPolyPatch::initGeometry()
{
    // Send over processorPolyPatch data
    processorPolyPatch::initGeometry();
}


void Foam::processorSplitCyclicPolyPatch::calcGeometry()
{
    // Receive and initialise processorPolyPatch data
    processorPolyPatch::calcGeometry();
Pout<< "processorSplitCyclicPolyPatch::calcGeometry()" << endl;
    // if (Pstream::parRun())
    // {

    //     // Where do we store the calculated transformation?
    //     // - on the processor patch?
    //     // - on the underlying cyclic patch?
    //     // - or do we not auto-calculate the transformation but
    //     //   have option of reading it.

    //     // Update underlying cyclic halves. Need to do both since only one
    //     // half might be present as a processorSplitCyclic.
    //     coupledPolyPatch& pp = const_cast<coupledPolyPatch&>(referPatch());
    //     pp.calcGeometry
    //     (
    //         *this,
    //         faceCentres(),
    //         faceAreas(),
    //         faceCellCentres(),
    //         neighbFaceCentres(),
    //         neighbFaceAreas(),
    //         neighbFaceCellCentres()
    //     );

    //     if (isA<splitCyclicPolyPatch>(pp))
    //     {
    //         const splitCyclicPolyPatch& cpp = refCast<const splitCyclicPolyPatch>(pp);
    //         const_cast<splitCyclicPolyPatch&>(cpp.neighbPatch()).calcGeometry
    //         (
    //             *this,
    //             neighbFaceCentres(),
    //             neighbFaceAreas(),
    //             neighbFaceCellCentres(),
    //             faceCentres(),
    //             faceAreas(),
    //             faceCellCentres()
    //         );
    //     }
    // }
}


void Foam::processorSplitCyclicPolyPatch::initMovePoints
(
    const pointField& p
)
{
    // Recalculate geometry
    initGeometry();
}


void Foam::processorSplitCyclicPolyPatch::movePoints
(
    const pointField&
)
{
    calcGeometry();
}


void Foam::processorSplitCyclicPolyPatch::initUpdateMesh()
{
    processorPolyPatch::initUpdateMesh();
}


void Foam::processorSplitCyclicPolyPatch::updateMesh()
{
     referPatchID_ = -1;
     processorPolyPatch::updateMesh();
}


void Foam::processorSplitCyclicPolyPatch::initOrder
(
    const primitivePatch& pp
) const
{
    // For now use the same algorithm as processorPolyPatch
    processorPolyPatch::initOrder(pp);
}


// Return new ordering. Ordering is -faceMap: for every face index
// the new face -rotation:for every new face the clockwise shift
// of the original face. Return false if nothing changes (faceMap
// is identity, rotation is 0)
bool Foam::processorSplitCyclicPolyPatch::order
(
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    // For now use the same algorithm as processorPolyPatch
    return processorPolyPatch::order( pp, faceMap, rotation);
}


void Foam::processorSplitCyclicPolyPatch::write(Ostream& os) const
{
    processorPolyPatch::write(os);
    os.writeKeyword("referPatch") << referPatchName_
        << token::END_STATEMENT << nl;
    if (tag_ != -1)
    {
        os.writeKeyword("tag") << tag_
            << token::END_STATEMENT << nl;
    }
}


// ************************************************************************* //

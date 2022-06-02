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

\*---------------------------------------------------------------------------*/

#include "splitCyclicPointPatch.H"
#include "pointBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "pointMesh.H"
#include "edgeList.H"
#include "transform.H"
#include "matchPoints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(splitCyclicPointPatch, 0);
    addToRunTimeSelectionTable
    (
        facePointPatch,
        splitCyclicPointPatch,
        polyPatch
    );
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //


void Foam::splitCyclicPointPatch::initGeometry()
{
    //polyPatch::initGeometry();
}


void Foam::splitCyclicPointPatch::calcGeometry()
{
    //polyPatch::calcGeometry();

    //calcTransforms();
}


/*void Foam::splitCyclicPointPatch::initGeometry(PstreamBuffers&)
{}


void Foam::splitCyclicPointPatch::calcGeometry(PstreamBuffers&)
{}
*/

void Foam::splitCyclicPointPatch::initMovePoints( const pointField&)
{}


void Foam::splitCyclicPointPatch::movePoints( )
{}


void Foam::splitCyclicPointPatch::initUpdateMesh( )
{
    facePointPatch::initUpdateMesh();
    splitCyclicPointPatch::initGeometry();
}


void Foam::splitCyclicPointPatch::updateMesh()
{
    facePointPatch::updateMesh();
    splitCyclicPointPatch::calcGeometry();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::splitCyclicPointPatch::splitCyclicPointPatch
(
    const polyPatch& patch,
    const pointBoundaryMesh& bm
)
:
    coupledFacePointPatch(patch, bm),
    splitCyclicPolyPatch_(refCast<const splitCyclicPolyPatch>(patch))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::splitCyclicPointPatch::~splitCyclicPointPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::edgeList& Foam::splitCyclicPointPatch::transformPairs() const
{
    return splitCyclicPolyPatch_.coupledPoints();
}


// ************************************************************************* //

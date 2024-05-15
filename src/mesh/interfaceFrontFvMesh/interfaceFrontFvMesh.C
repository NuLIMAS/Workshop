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

\*---------------------------------------------------------------------------*/

#include "interfaceFrontFvMesh.H"
#include "zeroGradientFvPatchFields.H"
#include "calculatedFvPatchFields.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "motionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "pointMesh.H"
#include "pointPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(interfaceFrontFvMesh, 0);

    addToRunTimeSelectionTable
    (
        topoChangerFvMesh,
        interfaceFrontFvMesh,
        IOobject
    );
}



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components

Foam::interfaceFrontFvMesh::interfaceFrontFvMesh(const IOobject& io)
:
    topoChangerFvMesh(io),
    motionDict_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).subDict(typeName + "Coeffs")
    ),
    wallFzID_(faceZones().findZoneID(motionDict_.lookup("wallFz"))),
    motionPtr_(motionSolver::New(*this))
    {
        Pout << "interfaceFrontFvMesh::interfaceFrontFvMesh" << endl;
        // if (openFzID_ < 0)
        // {
        //     FatalErrorIn("Foam::interfaceFrontFvMesh::interfaceFrontFvMesh(const IOobject& io)")
        //         << "open face zone does not exists." << nl
        //         << abort(FatalError);
        // }

        if (wallFzID_ < 0)
        {
            FatalErrorIn("Foam::interfaceFrontFvMesh::interfaceFrontFvMesh(const IOobject& io)")
                << "wall face zone does not exists." << nl
                << abort(FatalError);
        }

        topoChanger_.setSize(1);

        topoChanger_.set
        (
            0,
            new interfaceFront
            (
                "liq",
                0,
                topoChanger_,
                word(motionDict_.lookup("wallFz")),
                word(motionDict_.lookup("master")),
                word(motionDict_.lookup("slave")),
                scalarField(0, 0.0),
                true
            )
        );

        // Write mesh and modifiers
        topoChanger_.writeOpt() = IOobject::AUTO_WRITE;
        topoChanger_.write();
        write();

        faceZone& wall = faceZones()[wallFzID_];
        wall.resetAddressing(labelList(0), boolList(0, false));

    }


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::interfaceFrontFvMesh::~interfaceFrontFvMesh()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::interfaceFrontFvMesh::update()
{
    //Info << "interfaceFront 0" << endl;
    // Save old points
    pointField oldPointsNew = allPoints();
    //Info << "allPoints().size() " << allPoints().size() <<endl;
    tmp<pointField> newPoints = motionPtr_->newPoints();
    movePoints(newPoints());
    //movePoints(oldPointsNew); 

    const pointVectorField& pointMotionU = lookupObject<pointVectorField>("pointMotionU"); 
    //Info << "pointMotionU1 :" << pointMotionU << endl;

    // Get the interfaceFront object from topoChanger
    interfaceFront& attDet = dynamic_cast<interfaceFront&>(topoChanger_[0]);

    // Flag to check if the interface has changed
    //bool hasChanged = false;
    //Info << "interfaceFront 1" << endl;
    // Get the wall faceZone
    faceZone& wall = faceZones()[wallFzID_];

    pointField mappedOldPointsNew1;
    //Info << "State : " << attDet.attached() << endl;
    // Check if the interface is detached
    if ( !attDet.attached() )
    {
        // Set the attachment and change the mesh
        attDet.setAttach();
        autoPtr<mapPolyMesh> map = topoChanger_.changeMesh();
 
        motionPtr_->updateMesh(map());

        mappedOldPointsNew1.resize(allPoints().size());
        mappedOldPointsNew1.map(oldPointsNew, map->pointMap());

        movePoints(mappedOldPointsNew1);
        resetMotion();
        setV0();

        movePoints(map->preMotionPoints());
    }

    // Lookup whether the cell is  liqufied
    const volScalarField& iLiq = lookupObject<volScalarField>("Cl");

    // Get the neighbor and owner label lists
    const labelList& nei = neighbour();

    const labelList& own = owner();

    // Lists to store the new wall faces and their orientation
    DynamicList<label> newWall(nei.size());
    DynamicList<bool> newWallFlipMap(nei.size());
//Info << "newwall1:" << newWall<< endl;
    // Loop through neighbor faces
    forAll(nei, facei)
    {

    // Check if the neighbor is liquid and the owner is not, or vice versa
    if  (iLiq[nei[facei]] == 1 &&  iLiq[own[facei]] == 0)
    {
        newWall.append(facei);
        newWallFlipMap.append(true);
    }
    else if  (iLiq[nei[facei]] == 0 &&  iLiq[own[facei]] == 1)
    {
        newWall.append(facei);
        newWallFlipMap.append(false);
    }

    }
//Info << "newwall2:" << newWall<< endl;
    // Sort the new wall faces
    sort(newWall);

    // Reset addressing for the wall faceZone
    wall.resetAddressing(newWall, newWallFlipMap);

    pointField mappedOldPointsNew2;

    if (!mappedOldPointsNew1.size())
    {
        mappedOldPointsNew1 = allPoints();
    }
//Info << "pointMotionU2 :" << pointMotionU << endl;
    // Check if the interface is attached
    if ( attDet.attached() )
    {

        // Set the detachment and change the mesh
        attDet.setDetach();
        autoPtr<mapPolyMesh> map = topoChanger_.changeMesh();


        motionPtr_->updateMesh(map());
        //Info << "pointMotionU21 :" << pointMotionU << endl;


        //motionPtr_->updateMesh(map());
        //Info << "pointMotionU22 :" << pointMotionU << endl;    
        mappedOldPointsNew2.map(mappedOldPointsNew1, map->pointMap());
        movePoints(mappedOldPointsNew2);
        //Info << "pointMotionU23 :" << pointMotionU << endl;
        resetMotion();
        setV0();
 
        //movePoints(map->preMotionPoints());

    }
        //Info << "State : " << attDet.attached() << endl;
        //Info << "End of interfaceFrontFvMesh" << endl;
//Info << "pointMotionU3 :" << pointMotionU << endl;
    return false;
}


// ************************************************************************* //

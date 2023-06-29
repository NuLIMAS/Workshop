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

\*---------------------------------------------------------------------------*/

#include "simpleLiqFvMesh.H"
#include "zeroGradientFvPatchFields.H"
#include "calculatedFvPatchFields.H"
#include "surfaceFields.H"
#include "volFields.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(simpleLiqFvMesh, 0);

    addToRunTimeSelectionTable
    (
        topoChangerFvMesh,
        simpleLiqFvMesh,
        IOobject
    );
}



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components

Foam::simpleLiqFvMesh::simpleLiqFvMesh(const IOobject& io)
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
    //openFzID_(faceZones().findZoneID(motionDict_.lookup("openFz"))),
    wallFzID_(faceZones().findZoneID(motionDict_.lookup("wallFz")))
{
    if (openFzID_ < 0)
    {
        FatalErrorIn("Foam::simpleLiqFvMesh::simpleLiqFvMesh(const IOobject& io)")
            << "open face zone does not exists." << nl
            << abort(FatalError);
    }

    if (wallFzID_ < 0)
    {
        FatalErrorIn("Foam::simpleLiqFvMesh::simpleLiqFvMesh(const IOobject& io)")
            << "wall face zone does not exists." << nl
            << abort(FatalError);
    }

    topoChanger_.setSize(1);

    topoChanger_.set
    (
        0,
        new attachDetach
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
    //labelList newWall(wall);
    //sort(newWall);
    //wall.resetAddressing(newWall, Foam::boolList(newWall.size(), false));
    wall.resetAddressing(labelList(0), boolList(0, false));
        
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::simpleLiqFvMesh::~simpleLiqFvMesh()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::simpleLiqFvMesh::update()
{
    

    attachDetach& attDet = dynamic_cast<attachDetach&>(topoChanger_[0]);

    //const faceZone& open = faceZones()[openFzID_];
    //const faceZone& wall = faceZones()[wallFzID_];
    bool hasChanged = false;
    
    
    
    

    faceZone& wall = faceZones()[wallFzID_];


    if ( !attDet.attached() )
    {
        attDet.setAttach();
        autoPtr<mapPolyMesh> map = topoChanger_.changeMesh();
    }
    else
    {
    //    Info << "Already attached" << endl;
    }

    const label masterPatchIndex = attDet.masterPatchID().index();
    const volScalarField& iLiq = lookupObject<volScalarField>("Cl");



    //const scalarField& iLiq = liqueFlag.internalField();

    const labelList& nei = neighbour();

    const labelList& own = owner();

    
    DynamicList<label> newWall(nei.size());
    DynamicList<bool> newWallFlipMap(nei.size());
    forAll(nei, facei)
    {

    if  (iLiq[nei[facei]] == 1 &&  iLiq[own[facei]] == 0)
    {
        //Info << "Interface :" << facei << endl;
        newWall.append(facei);
        newWallFlipMap.append(true);
    }
    else if  (iLiq[nei[facei]] == 0 &&  iLiq[own[facei]] == 1)
    {
        //Info << "Interface1 :" << facei << endl;
        newWall.append(facei);
        newWallFlipMap.append(false);
    }
            
    }
    sort(newWall);

    wall.resetAddressing(newWall, newWallFlipMap);

    if ( attDet.attached() )
    {
    //    Info << "Detaching" << endl;
         attDet.setDetach();
         autoPtr<mapPolyMesh> map = topoChanger_.changeMesh();
    }
    else
    {
    //     Info << "Already detached" << endl;
    }




    {
        //faceZone& open = faceZones()[openFzID_];
        faceZone& wall = faceZones()[wallFzID_];

        Info << "End of simpleLiqFvMesh" << endl;
        //Info << open << endl;
        //Info << wall << endl;
    }
    return false;
}


// ************************************************************************* //

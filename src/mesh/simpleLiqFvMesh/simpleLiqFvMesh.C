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



// * * * * * * * * * * * * * Private Member Functions  * * *q * * * * * * * * //

/*
void Foam::simpleLiqFvMesh::addZonesAndModifiers()
{
    Info<< " Time = " << time().timeName() << endl
        << " Adding zones and modifiers to the mesh" << endl;
        

        

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
}

*/



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
    
    
    
    
    
    //const labelList& nei = mesh.neighbour();

    //const labelList& own = mesh.owner();

    
    
    

/*
    if ( !attDet.attached() )
    {
        const label masterPatchIndex = attDet.masterPatchID().index();
        const vectorField& Cfmp = Cf().boundaryField()[masterPatchIndex];
        forAll(Cfmp, i)
        {
            if(!liqFaceSelector.closed(Cfmp[i]))
            {
                Info << "a face has opened" << endl;
                hasChanged = true;
                break;
            }
        }

        if (!hasChanged)
        {
            forAll(open, oI)
            {
                if(liqFaceSelector.closed(Cf()[open[oI]]))
                {
                    Info << "a face has closed" << endl;
                    hasChanged = true;
                    break;
                }
            }
        }
    }
    else
    {
        hasChanged = true;
    }

    hasChanged = returnReduce(hasChanged, orOp<bool>());
*/

   // if (hasChanged)
   // {
    //    Info<< "Updating valve." << endl;

        //faceZone& open = faceZones()[openFzID_];
        faceZone& wall = faceZones()[wallFzID_];

        // Need to handle open faceZone by hand as it may overlap with DSMC
        //labelList openSave(open);
        //open.resetAddressing(labelList(), boolList());


        //Info << "Before setting wall zone" << endl;
        //Info << open << endl;
        //Info << wall << endl;

        if ( !attDet.attached() )
        {
      //      Info << "Attaching" << endl;

            attDet.setAttach();
            autoPtr<mapPolyMesh> map = topoChanger_.changeMesh();


          //  forAll(openSave, fzI)
          //  {
          //      openSave[fzI] = map->reverseFaceMap()[openSave[fzI]];
          //  }
            
        }
        else
        {
      //      Info << "Already attached" << endl;
        }

        //Info << "After Attaching" << endl;
        const label masterPatchIndex = attDet.masterPatchID().index();
       // Info << "# faces in open = " << returnReduce(open.size(), sumOp<label>())
      //  Info   << " # faces in wall = " << returnReduce(wall.size(), sumOp<label>())
      //      << " # faces in master patch = " << returnReduce(Cf().boundaryField()[masterPatchIndex].size(), sumOp<label>())
      //      << endl;
        //Info << open << endl;
        //Info << wall << endl;

/*
        DynamicList<label> newOpen(open.size() + wall.size());
        DynamicList<label> newWall(open.size() + wall.size());

        forAll(open, oI)
        {
            const label fI = open[oI];

            if(liqFaceSelector.closed(Cf()[fI]))
            {
                newWall.append(fI);
                
            }
            else
            {
                newOpen.append(fI);
            }
        }
        
        forAll(wall, wI)
        {
            const label fI = wall[wI];

           if(liqFaceSelector.closed(Cf()[fI]))
            {
                newWall.append(fI);
                 Info << "Wall" << fI << endl;
            }
            else
            {
                newOpen.append(fI);
               
            }
        }

        sort(newWall);
*/

//const fvPatchField<scalar>& iLiq =
//        patch().lookupPatchField<volScalarField, scalar>("liqueFlag");

const volScalarField& iLiq = lookupObject<volScalarField>("liqueFlag");



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
//Info << "newWall" << newWall << endl;


        wall.resetAddressing(newWall, newWallFlipMap);
        //wall.resetAddressing(labelList(0), boolList(0, false));
        //open.resetAddressing(newOpen, boolList(newOpen.size(), false));
        
            //wall.resetAddressing(labelList(0), boolList(0, false));
        

        //Info << "Before Detaching" << endl;
        //Info << open << endl;
        //Info << wall << endl;

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




        //forAll(newOpen, fzI)
        //{
        //    newOpen[fzI] = map->reverseFaceMap()[newOpen[fzI]];
        //}
        
        //open.resetAddressing(newOpen, boolList(newOpen.size(), false));

       //const label masterPatchIndex = attDet.masterPatchID().index();
      //  Info << "# faces in open = " << returnReduce(open.size(), sumOp<label>())
      //   Info   << " # faces in wall = " << returnReduce(wall.size(), sumOp<label>())
      //      << " # faces in master patch = " << returnReduce(Cf().boundaryField()[masterPatchIndex].size(), sumOp<label>())
      //      << endl;

   //correctPatches();

        // Trigger regeneration of some mesh addressing to avoid problems in parallel
        // Should be done in DSMC cloud, but there is no place to cover cell init
        // which happens before evolve.
 //      this->polyMesh::clearAddressing();
 //       this->tetBasePtIs();
 //       this->geometricD();


        //return true;
   // }
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

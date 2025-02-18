// Calculate the liqeuafaction criteria
liqueCriteria.internalField() = pE.internalField()/sigma0.internalField();
liqueCriteria.correctBoundaryConditions();

scalarField& ClI = Cl.internalField();
scalarField& CsI = Cs.internalField();

forAll (mesh.C(),cellI)
{
    if ( Cl[cellI] == 1 &&  Cs[cellI] == 0 )
    {
        ClI[cellI] = 1;
    }

    if ( liqueCriteria[cellI] >= 1.0 &&  Cs[cellI] == 0 )
    {
        ClI[cellI] = 1;
        //alpha[cellI] = n[cellI];
    }

    if (alpha[cellI]>= (1/(1+emin.value()))) //compacted  region
    {
        //Info  << " compacted" << endl;
        CsI[cellI] = 1;
        //E[cellI] = E2.value();
        ClI[cellI] = 0;
        alpha[cellI]=1/(1+emin.value());
        n[cellI]  =  emin.value()/(1+emin.value());
        relDensity[cellI] = relDensity2.value();

    }

}

Cl.correctBoundaryConditions();
Cs.correctBoundaryConditions();

// label patchID1 = mesh.boundaryMesh().findPatchID("sides");
// fvPatchScalarField& Clface1=Cl.boundaryField()[3];
// const labelList cl1=Clface1.patch().faceCells();
// Pout << "Cl " << Cl << endl;


// forAll(Clface1, facei)
// {
//   Pout << "cl :" <<  Cl[cl1[facei]] << endl;
// }
// If left and right faces doesn't liqueify together,
// switch face liquefied to non-liquefied


// Pout << "neignbor field" << Cl.boundaryField()[7] << endl;

scalarField Clneigh = Cl.boundaryField()[7].patchNeighbourField();
scalarField Clint = Cl.boundaryField()[7].patchInternalField();
const labelList cli=Cl.boundaryField()[7].patch().faceCells();


//  Pout <<  Cl.boundaryField()[7][1].owner() << endl;

forAll(Cl.boundaryField()[7],facei)
{

// //     Pout<< "values: "<< Clint[facei]  << " "
// //     << mesh.Cf().boundaryField()[7][facei] << " "
// //     <<  Clneigh[facei]
// //     << endl;
    if(Clneigh[facei] !=  Clint[facei])
    {
        Cl[cli[facei]] = 0;
    }
}


label patchID1 = mesh.boundaryMesh().findPatchID("left");
label patchID2 = mesh.boundaryMesh().findPatchID("right");

fvPatchScalarField& Clface1=Cl.boundaryField()[patchID1];
fvPatchScalarField& Clface2=Cl.boundaryField()[patchID2];
const labelList cl1=Clface1.patch().faceCells();
const labelList cl2=Clface2.patch().faceCells();
fvPatchScalarField& Csface1=Cs.boundaryField()[patchID1];
fvPatchScalarField& Csface2=Cs.boundaryField()[patchID2];
const labelList cs1=Csface1.patch().faceCells();
const labelList cs2=Csface2.patch().faceCells();
forAll(Clface1, facei)
{
    if ((Cl[cl1[facei]] != Cl[cl2[facei]] ) && ( (Cs[cs1[facei]] ==0)  || (Cs[cs2[facei]] ==0 ) ))
    {
        Cl[cl2[facei]] = 0;
        Cl[cl1[facei]] = 0;
    }
    else if  ((Cl[cl1[facei]] != Cl[cl2[facei]] ) && ( (Cs[cs1[facei]] != Cs[cs2[facei]])))
    {
        Cl[cl2[facei]] = 1;
        Cl[cl1[facei]] = 1;
    }

    //if (Cs[cs1[facei]] != Cs[cs2[facei]] )
    //{
    //    Cs[cs2[facei]] = 0;
    //    Cs[cs1[facei]] = 0;
    //}
}
forAll(Csface1, facei)
{
    if ((Cs[cs1[facei]] != Cs[cs2[facei]] ))
    {
        Cs[cs2[facei]] = 1;
        Cs[cs1[facei]] = 1;
    }
}


volScalarField waterHeight
(
    IOobject
    (
        "waterHeight",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    wh,
    zeroGradientFvPatchScalarField::typeName
);

// Initial effective stress
sigma0 = gammaD*(waterHeight-mesh.C().component(vector::Z))*((1.0+(2.0*k0))/3.0);

sigma0.correctBoundaryConditions();

liqueCriteria.internalField() = pE.internalField()/sigma0.internalField();
liqueCriteria.correctBoundaryConditions();

scalarField& ClI = Cl.internalField();

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
}

sigma0.correctBoundaryConditions();
Cl.correctBoundaryConditions();

// If left and right faces doesn't liqueify together, 
// switch face liquefied to non-liquefied    

label patchID1 = mesh.boundaryMesh().findPatchID("left");
label patchID2 = mesh.boundaryMesh().findPatchID("right");

fvPatchScalarField& Clface1=Cl.boundaryField()[patchID1];
fvPatchScalarField& Clface2=Cl.boundaryField()[patchID2];
const labelList cl1=Clface1.patch().faceCells();
const labelList cl2=Clface2.patch().faceCells();
forAll(Clface1, facei)
{
    if (Cl[cl1[facei]] != Cl[cl2[facei]] )
    {
        Cl[cl2[facei]] = 0;
        Cl[cl1[facei]] = 0;
    }
}


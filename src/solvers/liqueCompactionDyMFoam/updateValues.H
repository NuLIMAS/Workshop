
surfaceScalarField alphaF
(
    "alphaF",
    fvc::interpolate(alpha)
);
Vdj = V0 & mesh.Sf()/mesh.magSf()*pow(1-alphaF,q);


label boundaryID1 = mesh.boundaryMesh().findPatchID("sol");
label boundaryID2 = mesh.boundaryMesh().findPatchID("liq");
Vdj.boundaryField()[boundaryID1] == 0;
Vdj.boundaryField()[boundaryID2] == 0;
const fvPatchScalarField& patchUsol = Cl.boundaryField()[boundaryID1];
const fvPatchScalarField& patchUliq = Cl.boundaryField()[boundaryID2];
const labelList Clsol=patchUsol.patch().faceCells();
const labelList Clliq=patchUliq.patch().faceCells();

// Update viscosity and density

//visf = muf*(1+2.5*alphaF+10.05*alphaF*alphaF + 0.00273*exp(16.6*alphaF));
visf = muf + C* log(1-alphaF);

rho = ((1.0 - alpha) * rhoc + alpha * rhos) * Cl + (1-Cl) * rhoc;

//Update density, Cl at the boundary
forAll(patchUliq, i)
{
    Cl.boundaryField()[boundaryID2][i] = 1;
    Cl.boundaryField()[boundaryID1][i] = 0;
    rho.boundaryField()[boundaryID2][i] = 1000;
    rho.boundaryField()[boundaryID1][i] = 1000;
}

rhom = (1.0 - alpha) * rhoc + alpha * rhos;


labelHashSet soilCells;
forAll(Cl, cellI)
{

    if (Cl[cellI]==0 ) // non-liquefied region
    {
        soilCells.insert(cellI);

    }

    if (alpha[cellI]>= (1.0 - n2.value())) //compacted  region
    {
        //Info  << " compacted" << endl;
        soilCells.insert(cellI);
        Cs[cellI] = 1;
        Cl[cellI] = 0;
        n[cellI]  =  n2.value() ;
        alpha[cellI]=1.0 - n2.value();
        relDensity[cellI] = relDensity2.value();

    }

}

relDensity.correctBoundaryConditions();
E.correctBoundaryConditions();
nu.correctBoundaryConditions();
labelList soilLabels = soilCells.toc();
Cl.correctBoundaryConditions();
n.correctBoundaryConditions();
Cs.correctBoundaryConditions();
//rho.correctBoundaryConditions();
Vdj.correctBoundaryConditions();
mu = E/(2.0*(1.0 + nu));
lambda = nu*E/((1.0 + nu)*(1.0 - 2.0*nu));

Switch planeStress(materialProperties.lookup("planeStress"));

if (planeStress)
{
//    Info<< "Plane Stress\n" << endl;
    lambda = nu*E/((1.0 + nu)*(1.0 - nu));
}
else
{
//    Info<< "Plane Strain\n" << endl;
}
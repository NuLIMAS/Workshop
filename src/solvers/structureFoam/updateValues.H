
Switch stationary(materialProperties.lookup("stationary"));
label boundaryID3 = mesh.boundaryMesh().findPatchID("foundationSeabed");
//dimensionedScalar us(materialProperties.lookup("Us"));

//rope forces
//int tp = ceil( runTime.time().value() / T.value()) ;
scalar A =ropeForce.value();
scalar pi = 22/7;

const cellZone& cz1 = mesh.cellZones()["structureZones"];
forAll(cz1, czi)
{
    rhodg[cz1[czi]] = vector(0,0,-rhoag.value());
    E[cz1[czi]] = Ea.value();
    rhom[cz1[czi]] =rhoa.value();
    rho[cz1[czi]] =rhoa.value();
    Cf[cz1[czi]] =0;
    Cl[cz1[czi]] =0;
}

rho.boundaryField()[boundaryID3] = rhoa.value();
rhom.boundaryField()[boundaryID3] = rhoa.value();
const cellZone& cz2 = mesh.cellZones()["ropeZones1"];
forAll(cz2, czi)
{
    scalar force1 = A*Foam::sin(2.0*pi/T.value()*runTime.time().value() + 2*pi*mesh.C()[cz2[czi]][0]/wl.value()) ;
    Us[cz2[czi]] = vector(0,0,force1)/mesh.V()[czi];

}
const cellZone& cz3 = mesh.cellZones()["ropeZones2"];
forAll(cz3, czi)
{
    scalar force2 = A*Foam::sin(2.0*pi/T.value()*runTime.time().value() + 2*pi*mesh.C()[cz3[czi]][0]/wl.value()) ;
    Us[cz3[czi]] = vector(0,0,force2)/mesh.V()[czi];
}

E.correctBoundaryConditions();
rhom.correctBoundaryConditions();
Cf.correctBoundaryConditions();

// Update Shear modulus and Lames constant
mu = E/(2.0*(1.0 + nu));
lambda = nu*E/((1.0 + nu)*(1.0 - 2.0*nu));

Switch planeStress(materialProperties.lookup("planeStress"));

if (planeStress)
{
    lambda = nu*E/((1.0 + nu)*(1.0 - nu));
}

// Update Coefficient of consolidation
cv = (mu * k / gamma) * ((2-2*nu)/ ((1-2*nu)+(2-2*nu)*(n*mu/Kprime)));
cv.correctBoundaryConditions();

volScalarField Dp2
(
    "Dp2",
    (Kprime/n)
);

volScalarField Dp3
(
    "Dp3",
    (k/gamma)
);
Dp2.correctBoundaryConditions();
Dp3.correctBoundaryConditions();



forAll (mesh.C(),cellI)
{

    if ( liqueCriteria[cellI] >= 1.2 ||  Cf[cellI] == 0 ) //pressure buildp is turned off
    {
        Cf[cellI] = 0;
        cv[cellI] = 0;
        //Cl[cellI] = 1;
    }

}

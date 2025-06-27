
Switch stationary(materialProperties.lookup("stationary"));
label boundaryID3 = mesh.boundaryMesh().findPatchID("foundationSeabed");
//dimensionedScalar us(materialProperties.lookup("Us"));

//rope forces

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





surfaceScalarField alphaF
(
    "alphaF",
    fvc::interpolate(alpha)
);
//Info <<"alphaF" << alphaF << endl;
Vdj = V0 & mesh.Sf()/mesh.magSf()*pow(1-alphaF,q);


label boundaryID1 = mesh.boundaryMesh().findPatchID("sol");
label boundaryID2 = mesh.boundaryMesh().findPatchID("liq");
// Set drift-flux at the innterface to 0
Vdj.boundaryField()[boundaryID1] == 0; 
Vdj.boundaryField()[boundaryID2] == 0;
const fvPatchScalarField& patchUsol = Cl.boundaryField()[boundaryID1];
const fvPatchScalarField& patchUliq = Cl.boundaryField()[boundaryID2];
const labelList Clsol=patchUsol.patch().faceCells();
const labelList Clliq=patchUliq.patch().faceCells();

// Update viscosity and density

//visf = muf*(1+2.5*alphaF+10.05*alphaF*alphaF + 0.00273*exp(16.6*alphaF));
visf = muf + C* log(1-alphaF);
rhom = (1.0 - alpha) * rhoc + alpha * rhos;
rho = ((1.0 - alpha) * rhoc + alpha * rhos) * Cl + (1-Cl) * rhoc;

//Update density, Cl at the interface
forAll(patchUliq, i)
{
    Cl.boundaryField()[boundaryID2][i] = 1;
    Cl.boundaryField()[boundaryID1][i] = 0;
    rho.boundaryField()[boundaryID2][i] = 1000;
    rho.boundaryField()[boundaryID1][i] = 1000;
}

relDensity.correctBoundaryConditions();
E.correctBoundaryConditions();
nu.correctBoundaryConditions();

Cl.correctBoundaryConditions();
n.correctBoundaryConditions();
//Cs.correctBoundaryConditions();
Vdj.correctBoundaryConditions();


//From structureFoam

Switch stationary(materialProperties.lookup("stationary"));

//dimensionedScalar ms(materialProperties.lookup("ms"));

int tp = ceil( runTime.time().value() / T.value()) ;

const cellZone& cz1 = mesh.cellZones()["structureZones"];
forAll(cz1, czi)
{
    //rhodg[cz1[czi]] = vector(0,0,-rhoag.value());
    E[cz1[czi]] = Ea.value();
    rhom[cz1[czi]] =rhoa.value();
    rho[cz1[czi]] =rhoa.value();
    Cf[cz1[czi]] =0;
    Cl[cz1[czi]] =0;
    Cs[cz1[czi]] =1;
    alpha[cz1[czi]] =0.5;
    //E[cellI] = 1e10;
    //rhom[cellI] =1450;
    //Cf[cellI] =0;
    //Cl[cellI] =0;
}
const cellZone& cz2 = mesh.cellZones()["ropeZones1"];
forAll(cz2, czi)
{
    if (tp %2 == 0)
    {
        Us[cz2[czi]] = Ur1.value();
    }
    else
    {
        Us[cz2[czi]] = Ur2.value();
    }
}
const cellZone& cz3 = mesh.cellZones()["ropeZones2"];
forAll(cz3, czi)
{
    if (tp %2 == 0)
    {
        Us[cz3[czi]] = Ur2.value();
    }
    else
    {
        Us[cz3[czi]] = Ur1.value();
    }
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
//    Info<< "Plane Stress\n" << endl;
    lambda = nu*E/((1.0 + nu)*(1.0 - nu));
}

// Update Coefficient of consolidation
cv = (mu * k / gamma) * ((2-2*nu)/ ((1-2*nu)+(2-2*nu)*(n*mu/Kprime)));
cv.correctBoundaryConditions();

#if DEBUG
Info << "------ Material properties  --------" << endl;
Info << "Elastic modulus  :   max: " <<  max(E).value() <<  " Pa min: " << min(E).value() << " Pa" << endl;
Info << "Shear modulus    :   max: " <<  max(mu).value() <<  " Pa min: " << min(mu).value() << " Pa" << endl;
Info << "Lames constant   :   max: " <<  max(lambda).value() <<  " Pa min: " << min(lambda).value() << " Pa" << endl;
Info << "Porosity         :   max: " <<  max(n).value() <<  " min: " << min(n).value() << "" << endl;
Info << "Poisson ratio    :   max: " <<  max(nu).value() <<  " min: " << min(nu).value() << "" << endl;
Info << "Permeability     :   max: " <<  max(k).value() <<  " m/s min: " << min(k).value() << " m/s" << endl;
Info << "Density          :   max: " <<  max(rho).value() <<  " kg/m^3 min: " << min(rho).value() << " kg/m^3" << endl;

Info << "--------------------------------------------------" << endl;

Info << "---- Pressurebuildup properties ----" << endl;
Info << "relativeDesnity  :   max: " <<  max(relDensity).value() <<  " min: " << min(relDensity).value() << " " << endl;
Info << "Coefficient of lateral Earth pressure  : " <<  k0.value()  << endl;
Info << "Coefficient of consolidation  :   max: " <<  max(cv).value() <<  " m^2/s min: " << min(cv).value() << " m^2/s" << endl;
Info << "---------------------------------------------------" << endl;
#endif

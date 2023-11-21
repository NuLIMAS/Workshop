
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
rhom = (1.0 - alpha) * rhoc + alpha * rhos;
rho = ((1.0 - alpha) * rhoc + alpha * rhos) * Cl + (1-Cl) * rhoc;

//Update density, Cl at the boundary
forAll(patchUliq, i)
{
    Cl.boundaryField()[boundaryID2][i] = 1;
    Cl.boundaryField()[boundaryID1][i] = 0;
    rho.boundaryField()[boundaryID2][i] = 1000;
    rho.boundaryField()[boundaryID1][i] = 1000;
}


labelHashSet soilCells;
forAll(Cl, cellI)
{

    if (Cl[cellI]==0 ) // non-liquefied region
    {
        soilCells.insert(cellI);
    }
    else
    {
        E[cellI]=E2.value();
    }

    if (alpha[cellI]>= (1/(1+emin.value()))) //compacted  region
    {
        //Info  << " compacted" << endl;
        soilCells.insert(cellI);
        Cs[cellI] = 1;
        Cl[cellI] = 0;
        n[cellI]  =  emin.value()/(1+emin.value()) ;
        alpha[cellI]=1/(1+emin.value());
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

// Update Coefficient of consolidation
cv = (mu * k / gamma) * ((2-2*nu)/ ((1-2*nu)+(2-2*nu)*(n*mu/Kprime)));
cv.correctBoundaryConditions();

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

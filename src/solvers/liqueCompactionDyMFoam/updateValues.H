// Calculate the liqeuafction criteria
liqueCriteria.internalField() = pE.internalField()/sigma0.internalField();
liqueCriteria.correctBoundaryConditions();


// Update values in liquefied and non-liquefied regions
scalarField& ClI = Cl.internalField();
labelHashSet soilCells;
forAll (mesh.C(),cellI)
{
    if (( Cl[cellI] == 1 &&  Cs[cellI] == 0 ) || ( liqueCriteria[cellI] >= 1.0 &&  Cs[cellI] == 0 ))
    {
        ClI[cellI] = 1;
    }

    if (Cl[cellI]==0 ) // Non-liquefied region
    {
        soilCells.insert(cellI);
    }

    if (alpha[cellI]>= (1/(1+emin.value()))) //compacted  region
    {
        soilCells.insert(cellI);
        Cs[cellI] = 1;
        Cl[cellI] = 0;
        n[cellI]  = emin.value()/ (1+ emin.value());
        alpha[cellI] = 1/(1+emin.value());
        //relDensity[cellI] = relDensity2.value();
        //E[cellI]=E2.value();
    }

}
labelList soilLabels = soilCells.toc();
// If left and right faces doesn't liqueify together,
// switch face liquefied to non-liquefied

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
    if (Cl[cl1[facei]] != Cl[cl2[facei]] )
    {
        Cl[cl2[facei]] = 0;
        Cl[cl1[facei]] = 0;
    }
    if (Cs[cs1[facei]] != Cs[cs2[facei]] )
    {
        Cs[cs2[facei]] = 0;
        Cs[cs1[facei]] = 0;
    }
}


bool meshChanged = mesh.update();

// Update drift flux
surfaceScalarField alphaF
(
    "alphaF",
    fvc::interpolate(alpha)
);
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

// Set viscosities between the liquefied and non liquefied regions
visf.boundaryField()[boundaryID1] = 50;
visf.boundaryField()[boundaryID2] = 50;

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
Cs.correctBoundaryConditions();
Vdj.correctBoundaryConditions();
k.correctBoundaryConditions();

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


Info << "------ Material properties  --------" << endl;
Info << "Elastic modulus  :   max: " <<  max(E).value() <<  " Pa min: " << min(E).value() << " Pa" << endl;
Info << "Shear modulus    :   max: " <<  max(mu).value() <<  " Pa min: " << min(mu).value() << " Pa" << endl;
Info << "Lames constant   :   max: " <<  max(lambda).value() <<  " Pa min: " << min(lambda).value() << " Pa" << endl;
Info << "Porosity         :   max: " <<  max(n).value() <<  " min: " << min(n).value() << "" << endl;
Info << "Poisson ratio    :   max: " <<  max(nu).value() <<  " min: " << min(nu).value() << "" << endl;
Info << "Permeability     :   max: " <<  max(k).value() <<  " m/s min: " << min(k).value() << " m/s" << endl;
Info << "--------------------------------------------------" << endl;

Info << "---- Pressurebuildup properties ----" << endl;
Info << "relativeDesnity  :   max: " <<  max(relDensity).value() <<  " min: " << min(relDensity).value() << " " << endl;
Info << "Coefficient of lateral Earth pressure  : " <<  k0.value()  << endl;
Info << "Coefficient of consolidation  :   max: " <<  max(cv).value() <<  " m^2/s min: " << min(cv).value() << " m^2/s" << endl;
Info << "---------------------------------------------------" << endl;

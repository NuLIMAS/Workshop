Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volTensorField  gradU
(
    IOobject
    (
        "gradU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::grad(U)
);

//volTensorField gradU = fvc::grad(U);

Info<< "Reading field Vm\n" << endl;
volVectorField Vm
(
    IOobject
    (
        "Vm",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading pore pressure field p\n" << endl;

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading period averaged excess pore pressure  field pE\n" << endl;

volScalarField pE
(
    IOobject
    (
        "pE",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading inital effective stress\n" << endl;
volScalarField sigma0
(
    IOobject
    (
        "sigma0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    gammaD*(-mesh.C().component(vector::Z))*((1.0+(2.0*k0))/3.0),
    zeroGradientFvPatchScalarField::typeName
    //mesh
);
sigma0.correctBoundaryConditions();

/*Info<< "Reading inital effective stress\n" << endl;
volScalarField sigma1
(
    IOobject
    (
        "sigma1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    gammaD*(-mesh.C().component(vector::Z))*((1.0+(2.0*k0))/3.0),
    zeroGradientFvPatchScalarField::typeName
);
*/

volScalarField lambda("lambda", nu*E/((1.0 + nu)*(1.0 - 2.0*nu)));
volScalarField mu ("mu", E/(2.0*(1.0 + nu)));
lambda.correctBoundaryConditions();
mu.correctBoundaryConditions();

Info<< "Calculating stress field sigmaD\n" << endl;
volSymmTensorField sigmaD
(
    IOobject
    (
        "sigmaD",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mu*twoSymm(fvc::grad(U)) + lambda*(I*tr(fvc::grad(U)))
);

Info<< "Calculating explicit part of div(sigma) divSigmaExp\n" << endl;
volVectorField divSigmaExp
(
    IOobject
    (
        "divSigmaExp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(sigmaD)
);

if (compactNormalStress)
{
    divSigmaExp -= fvc::laplacian(2*mu + lambda, U, "laplacian(DD,U)");
}
else
{
    divSigmaExp -= fvc::div((2*mu + lambda)*fvc::grad(U), "div(sigmaD)");
}

Info<< "Initializing sigmaXZ field\n" << endl;
volScalarField tauXZ
(
    IOobject
    (
        "tauXZ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    sigmaD.component(symmTensor::XZ)
);

volScalarField tauXZPrime2Mean
(
    IOobject
    (
        "tauXZPrime2Mean",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0.0", p.dimensions()*p.dimensions(),0.0)
);

Info<< "Initializing strain field\n" << endl;
volSymmTensorField strain
(
    IOobject
    (
        "strain",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
     ),
     0.5*twoSymm(fvc::grad(U))

);

Info<< "Initializing volumetric strain field\n" << endl;
volScalarField volStrain
(
    IOobject
    (
        "volStrain",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    tr(strain)
);

Info<< "Initializing liquefaction criteria\n" << endl;
volScalarField liqueCriteria
(
    IOobject
    (
        "liqueCriteria",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField numCyc
(
    IOobject
    (
        "numCyc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

Info<< "Initialising field Vdj\n" << endl;

surfaceScalarField Vdj
(
    IOobject
    (
        "Vdj",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Initialising field cv\n" << endl;
volScalarField cv
(
    IOobject
    (
        "cv",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0.0", dimLength*dimLength/dimTime,0.0)
);

Info<< "Initialising field f\n" << endl;
volScalarField f
(
    IOobject
    (
        "f",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0.0", dimPressure/dimTime,0.0)
);

// For piso loop
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p,piso.dict(), pRefCell, pRefValue);
mesh.schemesDict().setFluxRequired(p.name());

Info<< "Initialising field phi, phiAlpha\n" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    linearInterpolate(rho*Vm) & mesh.Sf()
);

surfaceScalarField phiAlpha
(
    IOobject
    (
        "phiAlpha",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0.0", dimVolume/dimTime,0.0)
);

volScalarField contErr
(
    IOobject
    (
        "contErr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

// Set initial  alpha and water height
forAll (mesh.C(),cellI)
{
    if (mesh.C().internalField()[cellI].component(vector::Z)> wh.value())
    {
        Cl[cellI] = 1;
        //E[cellI] = 1e10;
        alpha[cellI] = 0;
        Cf[cellI] = 0;
        //rhom[cellI] =1395;
    }
    //else
    //{
    //    Cl[cellI] = 1;
    //}

}



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

volVectorField Us
(
    IOobject
    (
        "Us",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimForce/dimVolume,vector (0.0, 0.0, 0)),
    zeroGradientFvPatchScalarField::typeName
);


// Initial effective stress
//sigma0 = gammaD*(waterHeight-mesh.C().component(vector::Z))*((1.0+(2.0*k0))/3.0);
//sigma0 = gammaD*(-mesh.C().component(vector::Z))*((1.0+(2.0*k0))/3.0);
volScalarField sigmaA
(
    IOobject
    (
        "sigmaA",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimForce/dimArea, 0.0),
    zeroGradientFvPatchScalarField::typeName

);
volScalarField dist
(
    IOobject
    (
        "dist",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimLength, 0),
    zeroGradientFvPatchScalarField::typeName
);
volScalarField fd
(
    IOobject
    (
        "fd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimLength, fw.value()/2),
    zeroGradientFvPatchScalarField::typeName
);

if(loading)
{
    label patchID = mesh.boundaryMesh().findPatchID("seabedFoundation");
    const vectorField& faceCentres = mesh.Cf().boundaryField()[patchID];
    forAll(mesh.C(),cellI)
    {
        scalar mink= GREAT;
        forAll(faceCentres, facei)
        {
            if ( mag(mesh.C()[cellI] - faceCentres[facei]) < mink)
            {
                mink =mag(mesh.C()[cellI] - faceCentres[facei]) ;
            }
        }
        dist[cellI] = mink;
    }

sigmaA.internalField() =QbyBL.value()* (1 -
    ((pow(dist.internalField(),3.0))/
    (pow(pow(fd.internalField(),2)+ pow(dist.internalField(),2),1.5))));
}

sigma0 = sigma0 + sigmaA * (1.0+(2.0*k0))/3.0;



#if DEBUG
Info << "Initial parameters used :" << endl;
Info << "----- Wave properties -----" << endl;
Info << "Wave height  : " <<  wahe.value() << " m" << endl;
Info << "Wave length  : " <<  wl.value() << " m" << endl;
Info << "Water depth  : " <<  wd.value() << " m" << endl;
Info << "Wave period  : " <<  T.value() << " s" << endl;
Info << "---------------------------------------" << endl;

Info << "-------- Material properties  --------" << endl;
Info << "Elastic modulus  :   max: " <<  max(E).value() <<  " Pa min: " << min(E).value() << " Pa" << endl;
Info << "Porosity         :   max: " <<  max(n).value() <<  " min: " << min(n).value() << "" << endl;
Info << "Poisson ratio    :   max: " <<  max(nu).value() <<  " min: " << min(nu).value() << "" << endl;
Info << "Permeability     :   max: " <<  max(k).value() <<  " m/s min: " << min(k).value() << " m/s" << endl;
Info << "-----------------------------------------" << endl;

Info << "-------- Pressurebuildup properties --------" << endl;
Info << "relativeDesnity  :   max: " <<  max(relDensity).value() <<  " min: " << min(relDensity).value() << " " << endl;
Info << "Max void ratio   : " <<  emax.value()  << endl;
Info << "Min void ratio   : " <<  emin.value()  << endl;
Info << "Coefficient of lateral Earth pressure  : " <<  k0.value()  << endl;
Info << "-----------------------------------------" << endl;

Info << "-------- Drift flux model properties --------" << endl;
Info << "Settling velocity  : " <<  V0.value()  << endl;
Info << "Exponent q (Vdj= V0*(1-alpha)^q)  : " <<  q.value()  << endl;
Info << "-----------------------------------------" << endl;

#endif

Info<< "Reading material properties\n" << endl;

IOdictionary materialProperties
(
    IOobject
    (
        "materialProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);


IOdictionary waveProperties
(
    IOobject
    (
        "waveProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);


Info<< "Initializing field E\n" << endl;
volScalarField E
(
    IOobject
    (
        "E",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(materialProperties.lookup("E")),
    zeroGradientFvPatchScalarField::typeName
);

Info<< "Initializing field k\n" << endl;
volScalarField k
(
    IOobject
    (
        "k",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(materialProperties.lookup("k")),
    zeroGradientFvPatchScalarField::typeName
);

Info<< "Initializing field nu\n" << endl;
volScalarField nu
(
    IOobject
    (
        "nu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(materialProperties.lookup("nu")),
    zeroGradientFvPatchScalarField::typeName
);

Info<< "Initialising field n\n" << endl;
volScalarField n
(
    IOobject
    (
        "n",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(materialProperties.lookup("n")),
    zeroGradientFvPatchScalarField::typeName
);
dimensionedScalar Ea(materialProperties.lookup("Ea"));
dimensionedScalar K(materialProperties.lookup("K"));
dimensionedScalar Sr(materialProperties.lookup("Sr"));
dimensionedScalar p0(materialProperties.lookup("p0"));
dimensionedScalar gamma(materialProperties.lookup("gamma"));
dimensionedScalar rhoc(materialProperties.lookup("rhof"));
dimensionedScalar rhos(materialProperties.lookup("rhos"));

//dimensionedVector g(materialProperties.lookup("g"));
dimensionedVector Ur1(materialProperties.lookup("Ur1"));
dimensionedVector Ur2(materialProperties.lookup("Ur2"));
dimensionedScalar rhoag(materialProperties.lookup("rhodg"));

dimensionedScalar wh(materialProperties.lookup("wh"));
dimensionedScalar k0(materialProperties.lookup("k0"));
dimensionedScalar gammaDash(materialProperties.lookup("gammaD"));
dimensionedScalar T("T",dimensionSet(0, 0, 1, 0, 0), waveProperties.lookup("period") );
dimensionedScalar wahe("wahe",dimensionSet(0, 1, 0, 0, 0), waveProperties.lookup("waveHeight"));
dimensionedScalar wl("wl",dimensionSet(0, 1, 0, 0, 0), waveProperties.lookup("waveLength"));
dimensionedScalar wd("wd",dimensionSet(0, 1, 0, 0, 0), waveProperties.lookup("waterDepth"));
Switch st(materialProperties.lookup("steady"));

Switch loading(materialProperties.lookup("loading"));
dimensionedScalar QbyBL(materialProperties.lookup("QbyBL"));
dimensionedScalar emax(materialProperties.lookup("emax"));
dimensionedScalar emin(materialProperties.lookup("emin"));
dimensionedScalar rhoa(materialProperties.lookup("rhoa"));
dimensionedScalar fw(materialProperties.lookup("fw"));

dimensionedScalar ropeForce(materialProperties.lookup("ropeForce"));
Info<< "Initialising field relDensity\n" << endl;
volScalarField relDensity
(
    IOobject
    (
        "relDensity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
     ),
    mesh,
    dimensionedScalar (materialProperties.lookup("relDensity")),
    zeroGradientFvPatchScalarField::typeName
);

Info<< "Initialising field gammaD\n" << endl;
volScalarField gammaD
(
    IOobject
    (
        "gammaD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
     ),
    mesh,
    dimensionedScalar (materialProperties.lookup("gammaD")),
    zeroGradientFvPatchScalarField::typeName
);


volVectorField rhodg
(
    IOobject
    (
        "rhodg",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimForce/dimVolume, vector(0, 0, 0)),
    zeroGradientFvPatchScalarField::typeName
);
rhodg.correctBoundaryConditions();

dimensionedScalar Kprime = 1.0/(1.0/K + (1.0-Sr)/p0);

Info<< "Calculating Lame's coefficients\n" << endl;


Info<< "Initializing field mu\n" << endl;
volScalarField mu
(
    IOobject
    (
         "mu",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet (1,-1,-2, 0, 0, 0, 0),
    zeroGradientFvPatchScalarField::typeName
);

Info<< "Initializing field lambda\n" << endl;
volScalarField lambda
(
    IOobject
    (
        "lambda",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet (1,-1,-2, 0, 0, 0, 0),
    zeroGradientFvPatchScalarField::typeName

);


Switch planeStress(materialProperties.lookup("planeStress"));

if (planeStress)
{
    Info<< "Plane Stress\n" << endl;
    lambda = nu*E/((1.0 + nu)*(1.0 - nu));
}
else
{
    Info<< "Plane Strain\n" << endl;
}

mu.correctBoundaryConditions();
lambda.correctBoundaryConditions();
k.correctBoundaryConditions();
nu.correctBoundaryConditions();

Info<< "Initialising field Cl\n" << endl;
volScalarField Cl
(
    IOobject
    (
        "Cl",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Cl",dimless,0),
    zeroGradientFvPatchScalarField::typeName
);


volScalarField Cf
(
    IOobject
    (
        "Cf",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Cf",dimless,1),
    zeroGradientFvPatchScalarField::typeName
);


Info<< "Initializing field rho, rhof, rhom\n" << endl;
volScalarField rhof
(
    IOobject
    (
        "rhof",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(materialProperties.lookup("rhof")),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    (n* rhoc + (1-n) *  rhos)*Cl + (1-Cl) * rhoc
);

volScalarField rhom
(
    IOobject
    (
        "rhom",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    (n * rhoc + (1-n) *  rhos)
);

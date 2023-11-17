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
        IOobject::NO_READ,
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
    dimensionedScalar(materialProperties.lookup("n1")),
    zeroGradientFvPatchScalarField::typeName
);


dimensionedScalar K(materialProperties.lookup("K"));
dimensionedScalar Sr(materialProperties.lookup("Sr"));
dimensionedScalar p0(materialProperties.lookup("p0"));
dimensionedScalar gamma(materialProperties.lookup("gamma"));
dimensionedScalar muf(materialProperties.lookup("vis"));
dimensionedScalar C(materialProperties.lookup("C"));
dimensionedScalar rhoc(materialProperties.lookup("rhof"));
dimensionedScalar rhos(materialProperties.lookup("rhos"));
dimensionedVector V0(materialProperties.lookup("V0"));
dimensionedVector g(materialProperties.lookup("g"));
dimensionedScalar q(materialProperties.lookup("q"));
dimensionedScalar wh(materialProperties.lookup("wh"));
dimensionedScalar k0(materialProperties.lookup("k0"));
dimensionedScalar gammaDash(materialProperties.lookup("gammaD"));
dimensionedScalar T("T",dimensionSet(0, 0, 1, 0, 0), waveProperties.lookup("period") );
dimensionedScalar wahe("wahe",dimensionSet(0, 1, 0, 0, 0), waveProperties.lookup("waveHeight"));
dimensionedScalar wl("wl",dimensionSet(0, 1, 0, 0, 0), waveProperties.lookup("waveLength"));
dimensionedScalar wd("wd",dimensionSet(0, 1, 0, 0, 0), waveProperties.lookup("waterDepth"));
//dimensionedScalar E2(materialProperties.lookup("E2"));
dimensionedScalar relDensity2(materialProperties.lookup("relDensity2"));
dimensionedScalar emax(materialProperties.lookup("emax"));
dimensionedScalar emin(materialProperties.lookup("emin"));
//volScalarField e = n /(1.0 -n);
Switch st(materialProperties.lookup("steady"));
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


mu.correctBoundaryConditions();
lambda.correctBoundaryConditions();
k.correctBoundaryConditions();
nu.correctBoundaryConditions();


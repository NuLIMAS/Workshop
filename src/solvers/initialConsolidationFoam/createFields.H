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
        IOobject::AUTO_WRITE
    ),
    fvc::grad(U)
);

Info<< "Calculating stress field sigmaD\n" << endl;
volSymmTensorField sigmaD
(
    IOobject
    (
        "sigmaD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
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
        IOobject::AUTO_WRITE
    ),
    fvc::div(sigmaD)
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
        IOobject::NO_WRITE
     ),
     0.5*twoSymm(fvc::grad(U))
);

volSymmTensorField currentStrain
(
    IOobject
    (
        "currentStrain",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
     ),
     strain
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
        IOobject::AUTO_WRITE
    ),
    tr(strain)
);

volVectorField rhog
(
    IOobject
    (
        "rhog",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimForce/dimVolume, vector(0, 0, -gammaD.value())),
    zeroGradientFvPatchScalarField::typeName
);
rhog.correctBoundaryConditions();

volScalarField effStress
(
    IOobject
    (
        "effStress",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (gammaD* -mesh.C().component(vector::Z))*
    (1 + 2*k0)/3
);

volScalarField sigmaA
(
    IOobject
    (
        "sigmaA",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimForce/dimArea, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

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
    mesh,
    dimensionedScalar("zero", dimForce/dimArea, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

divSigmaExp -= fvc::laplacian(2*mu + lambda, U, "laplacian(DD,D)");


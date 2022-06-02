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


dimensionedScalar p0(materialProperties.lookup("p0"));
dimensionedScalar gamma(materialProperties.lookup("gamma"));

Info<< "Initializing field E\n" << endl;
volScalarField E
(
    IOobject
    (
        "E",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
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
        IOobject::AUTO_WRITE
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
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(materialProperties.lookup("nu")),
    zeroGradientFvPatchScalarField::typeName
);

Info<< "Initializing field K\n" << endl;
volScalarField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(materialProperties.lookup("K")),
    zeroGradientFvPatchScalarField::typeName
);

Info<< "Initializing field n\n" << endl;
volScalarField n
(
    IOobject
    (
        "n",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(materialProperties.lookup("n")),
    zeroGradientFvPatchScalarField::typeName
);

Info<< "Initializing field Sr\n" << endl;
volScalarField Sr
(
    IOobject
    (
        "Sr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(materialProperties.lookup("Sr")),
    zeroGradientFvPatchScalarField::typeName
);



E.correctBoundaryConditions();
n.correctBoundaryConditions();
k.correctBoundaryConditions();
nu.correctBoundaryConditions();
K.correctBoundaryConditions();
Sr.correctBoundaryConditions();



volScalarField Kprime = 1.0/(1.0/K + (1.0-Sr)/p0);
volScalarField Kdash = 1.0/Kprime;

Info<< "Calculating Lame's coefficients\n" << endl;

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
    E/(2.0*(1.0 + nu))
);

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
    nu*E/((1.0 + nu)*(1.0 - 2.0*nu))
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





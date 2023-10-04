Info<< "Reading material properties\n" << endl;

IOdictionary materialProperties
(
    IOobject
    (
        "materialProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

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

dimensionedScalar rhos(materialProperties.lookup("rhos"));
dimensionedScalar rhof(materialProperties.lookup("rhof"));
dimensionedScalar g(materialProperties.lookup("g"));
dimensionedScalar k0(materialProperties.lookup("k0"));
dimensionedScalar n(materialProperties.lookup("n"));
dimensionedScalar gammaD("gammaD", (1-n)*(rhos-rhof)*g );
Switch planeStress(materialProperties.lookup("planeStress"));

Info<< "Calculating Shear modulus and Lame's coefficients\n" << endl;

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
    E/(2.0*(1.0 + nu)),
    zeroGradientFvPatchScalarField::typeName
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
    nu*E/((1.0 + nu)*(1.0 - 2.0*nu)),
    zeroGradientFvPatchScalarField::typeName
);

if (planeStress)
{
    Info<< "Plane Stress\n" << endl;
    lambda = nu*E/((1.0 + nu)*(1.0 - nu));
}
else
{
    Info<< "Plane Strain\n" << endl;
}

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

//dimensionedScalar E(materialProperties.lookup("E"));
//dimensionedScalar nu(materialProperties.lookup("nu"));
//dimensionedScalar k(materialProperties.lookup("k"));
dimensionedScalar K(materialProperties.lookup("K"));
dimensionedScalar Sr(materialProperties.lookup("Sr"));
dimensionedScalar p0(materialProperties.lookup("p0"));
dimensionedScalar gamma(materialProperties.lookup("gamma"));
//dimensionedScalar n(materialProperties.lookup("n"));

dimensionedScalar k0(materialProperties.lookup("k0"));
dimensionedScalar emax(materialProperties.lookup("emax"));
dimensionedScalar emin(materialProperties.lookup("emin"));
dimensionedScalar gammaD(materialProperties.lookup("gammaD"));
dimensionedScalar g("g", dimLength/dimTime/dimTime, 9.81);

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

k.correctBoundaryConditions();
nu.correctBoundaryConditions();

Info<< "Initializing field n\n" << endl;
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



dimensionedScalar Kprime = 1.0/(1.0/K + (1.0-Sr)/p0);



Info<< "Initialising field relDensity\n" << endl;
volScalarField relDensity
(
    IOobject
    (
        "relDensity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
     ),
    mesh,
    dimensionedScalar (materialProperties.lookup("relDensity")),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField Cl
(
    IOobject
    (
        "Cl",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Cl",dimless,0),
    zeroGradientFvPatchScalarField::typeName
);


volScalarField liqueFlag
(
    IOobject
    (
        "liqueFlag",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("liqueFlag",dimless,0)
);


Switch stones(materialProperties.lookup("stones"));
if(stones)
{
    dimensionedScalar Es(materialProperties.lookup("Es"));
    // dimensionedScalar nus(materialProperties.lookup("nus"));
    dimensionedScalar ns(materialProperties.lookup("ns"));
    dimensionedScalar relDensitys(materialProperties.lookup("relDensitys"));
    dimensionedScalar ks(materialProperties.lookup("ks"));

    const cellZone& cz1 = mesh.cellZones()["stone"];

    forAll(cz1, czi)
    {
        E[cz1[czi]] = Es.value();
        // nu[cz1[czi]] = nus.value();
        n[cz1[czi]] = ns.value();
        relDensity[cz1[czi]] = relDensitys.value();
        k[cz1[czi]] = ks.value();

    }

}

E.correctBoundaryConditions();
nu.correctBoundaryConditions();
k.correctBoundaryConditions();
relDensity.correctBoundaryConditions();

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
         IOobject::NO_WRITE
    ),
    E/(2.0*(1.0 + nu))
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
        IOobject::NO_WRITE
    ),
    nu*E/((1.0 + nu)*(1.0 - 2.0*nu))
);


mu.correctBoundaryConditions();
lambda.correctBoundaryConditions();

volScalarField Dp
(
    IOobject
    (
        "Dp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (k/gamma*Kprime/n),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField Dp2
(
    IOobject
    (
        "Dp2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (Kprime/n),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField Dp3
(
    IOobject
    (
        "Dp3",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (k/gamma),
    zeroGradientFvPatchScalarField::typeName
);


Info<< "Reading wave properties\n" << endl;

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


dimensionedScalar T("T",dimensionSet(0, 0, 1, 0, 0), waveProperties.lookup("period") );

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
dimensionedScalar g(materialProperties.lookup("g"));
dimensionedScalar vis(materialProperties.lookup("vis"));
dimensionedScalar rhof(materialProperties.lookup("rhof"));
dimensionedScalar rhos(materialProperties.lookup("rhos"));

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
        IOobject::NO_WRITE
    ),    
    mesh,
    dimensionedScalar(materialProperties.lookup("nu")),
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
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(materialProperties.lookup("n")),
    zeroGradientFvPatchScalarField::typeName
);



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

volScalarField Dp4
(
    IOobject
    (
        "Dp4",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (k/g),
    zeroGradientFvPatchScalarField::typeName
);


Info<< "Initializing field Dp5\n" << endl;
volScalarField Dp5
(
    IOobject
    (
        "Dp5",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("one", dimless, 1.0),
    zeroGradientFvPatchScalarField::typeName
); 



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



dimensionedScalar relDensity(materialProperties.lookup("relDensity"));
dimensionedScalar T("T",dimensionSet(0, 0, 1, 0, 0), waveProperties.lookup("period") );
dimensionedScalar gammaD(materialProperties.lookup("gammaD"));
dimensionedScalar alpha (materialProperties.lookup("alpha"));
dimensionedScalar beta (materialProperties.lookup("beta"));
dimensionedScalar k0 (materialProperties.lookup("k0"));
//dimensionedScalar alpha = 0.34 * relDensity + 0.084;
//dimensionedScalar beta = 0.37 * relDensity - 0.46;




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
    gammaD*(max(mesh.C().component(vector::Z)) -    mesh.C().component(vector::Z))*((1.0+(2.0*k0))/3.0)
);

volScalarField tauAmp = sqrt(2.0* tauRMS);

liqueCriteria.internalField() = pE.internalField()/sigma0.internalField();
liqueCriteria.correctBoundaryConditions();


volScalarField numCyc 
(
    IOobject
    (
        "numCyc",
        runTime.timeName(),
        mesh, 
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);


numCyc.internalField() =
    Foam::pow
    (
        (1.0/alpha.value())
       *mag(tauAmp.internalField()/sigma0.internalField()) + SMALL, (1.0/beta.value())
    );

   
numCyc.correctBoundaryConditions();

Info<< "Calculate source term\n" << endl;   
        
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
    sigma0/(numCyc * T)
);
    
f.correctBoundaryConditions();
runTime.write();  
   



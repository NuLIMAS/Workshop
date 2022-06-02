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

dimensionedScalar E(materialProperties.lookup("E"));
dimensionedScalar nu(materialProperties.lookup("nu"));
dimensionedScalar k(materialProperties.lookup("k"));
dimensionedScalar K(materialProperties.lookup("K"));
dimensionedScalar Sr(materialProperties.lookup("Sr"));
dimensionedScalar p0(materialProperties.lookup("p0"));
dimensionedScalar gamma(materialProperties.lookup("gamma"));
dimensionedScalar n(materialProperties.lookup("n"));

dimensionedScalar Kprime = 1.0/(1.0/K + (1.0-Sr)/p0);

Info<< "Calculating Lame's coefficients\n" << endl;

dimensionedScalar mu = E/(2.0*(1.0 + nu));
dimensionedScalar lambda = nu*E/((1.0 + nu)*(1.0 - 2.0*nu));

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

Info<< "mu = " << mu.value() << " Pa\n";
Info<< "lambda = " << lambda.value() << " Pa\n";

// Fake some fields
volScalarField muField
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    mu
);

volScalarField lambdaField
(
    IOobject
    (
        "lambda",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    lambda
);


dimensionedScalar Dp
(
    "Dp",
    dimensionSet(0, 2, -1 , 0, 0),
    (k/gamma*Kprime/n).value()
);

dimensionedScalar Dp2
(
    "Dp2",
    dimensionSet(1, -1, -2 , 0, 0),
    (Kprime/n).value()
);

dimensionedScalar Dp3
(
    "Dp3",
    dimensionSet(-1, 3, 1 , 0, 0),
    (k/gamma).value()
);
dimensionedScalar Dp4
(
    "Dp4",
    dimensionSet(0, 0, 1 , 0, 0),
    (k/9.81).value()
);

Info<< "Dp = " << Dp.value() << " m^2/s \n";
Info<< "Dp2 = " << Dp2.value() << " kg/m/s^2\n";
Info<< "Dp3 = " << Dp3.value() << "m^3/kg\n";
Info<< "Dp4 = " << Dp3.value() << "s\n";
Info << "Fluid bulk mod " << Kprime << endl;




Info<< "Initializing field Cs\n" << endl;
volScalarField Cs
(
    IOobject
    (
        "Cs",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("one", dimless, 1.0)
);

Info<< "creating liquefaction flag field liqueFlag\n" << endl;

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
    dimensionedScalar("liqueFlag",dimless,1.0)
); 



/*

const cellZone& cz1 = mesh.cellZones()["leftZones"];
forAll(cz1, czi)
{
//Info <<"I am here : left"<< endl;
Cs[cz1[czi]] = 0;

}


const cellZone& cz2 = mesh.cellZones()["rightZones"];
forAll(cz2, czi)
{
//Info <<"I am here : right"<< endl;
Cs[cz2[czi]] = 0;

}

const cellZone& cz3 = mesh.cellZones()["rockerZones"];
forAll(cz3, czi)
{

Cs[cz3[czi]] = 0;

}
*/  
  


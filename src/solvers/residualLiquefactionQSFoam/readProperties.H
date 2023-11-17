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

dimensionedScalar E(materialProperties.lookup("E"));
dimensionedScalar nu(materialProperties.lookup("nu"));
dimensionedScalar k(materialProperties.lookup("k"));
dimensionedScalar K(materialProperties.lookup("K"));
dimensionedScalar Sr(materialProperties.lookup("Sr"));
dimensionedScalar p0(materialProperties.lookup("p0"));
dimensionedScalar gamma(materialProperties.lookup("gamma"));
dimensionedScalar n(materialProperties.lookup("n"));
dimensionedScalar e = n/(1.0 -n);
dimensionedScalar gammaD(materialProperties.lookup("gammaD"));
dimensionedScalar k0(materialProperties.lookup("k0"));
dimensionedScalar relDensity(materialProperties.lookup("relDensity"));

dimensionedScalar T("T",dimensionSet(0, 0, 1, 0, 0), waveProperties.lookup("period") );
dimensionedScalar wahe("wahe",dimensionSet(0, 1, 0, 0, 0), waveProperties.lookup("waveHeight"));
dimensionedScalar wl("wl",dimensionSet(0, 1, 0, 0, 0), waveProperties.lookup("waveLength"));
dimensionedScalar wd("wd",dimensionSet(0, 1, 0, 0, 0), waveProperties.lookup("waterDepth"));

dimensionedScalar Kprime = 1.0/(1.0/K + (1.0-Sr)/p0);

// Calculate empirical coefficients alpha and beta
dimensionedScalar alpha = 0.34 * relDensity + 0.084;
dimensionedScalar beta = 0.37 * relDensity - 0.46;

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

Info << "Parameters used :" << endl;
Info << "----- Wave properties -----" << endl;
Info << "Wave height  : " <<  wahe.value() << " m" << endl;
Info << "Wave length  : " <<  wl.value() << " m" << endl;
Info << "Water depth  : " <<  wd.value() << " m" << endl;
Info << "Wave period  : " <<  T.value() << " s" << endl;
Info << "---------------------------------------" << endl;

Info << "------ Material properties  --------" << endl;
Info << "Elastic modulus  :   max: " <<  E.value() <<  " Pa "<< endl;
Info << "Shear modulus    :   max: " <<  mu.value() <<  " Pa "<< endl;
Info << "Lames constant   :   max: " <<  lambda.value() <<  " Pa "<< endl;
Info << "Porosity         :   max: " <<  n.value() << endl;
Info << "Poisson ratio    :   max: " <<  nu.value() << endl;
Info << "Permeability     :   max: " <<  k.value() <<  " m/s "<< endl;
Info << "--------------------------------------------------" << endl;
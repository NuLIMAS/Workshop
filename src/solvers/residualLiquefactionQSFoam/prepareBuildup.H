

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




dimensionedScalar e = n/(1.0 -n);
dimensionedScalar k0(materialProperties.lookup("k0"));
dimensionedScalar relDensity(materialProperties.lookup("relDensity"));
//dimensionedScalar alpha(materialProperties.lookup("alpha"));
//dimensionedScalar beta(materialProperties.lookup("beta"));
dimensionedScalar alpha = 0.34 * relDensity + 0.084;
dimensionedScalar beta = 0.37 * relDensity - 0.46;

dimensionedScalar T("T",dimensionSet(0, 0, 1, 0, 0), waveProperties.lookup("period") );

volScalarField tauXY("tauXY", sigmaD.component(symmTensor::XY));
tauXZ = sigmaD.component(symmTensor::XZ);
volScalarField tauYZ("tauYZ", sigmaD.component(symmTensor::YZ));

dimensionedScalar gammaD(materialProperties.lookup("gammaD"));
Info<< "runTime :" << runTime.time().value() <<endl;

volScalarField txz = SMALL*sigmaD.component(symmTensor::XZ);
if ( runTime.time().value() > T.value())
{
     txz = sqrt(mesh.lookupObject<volScalarField>("tauXZPrime2Mean")*2.0);
}

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

volScalarField tauEq
(
    IOobject
    (
        "tauEq",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    sqrt((tauXY*tauXY)+(tauXZ*tauXZ)+(tauYZ*tauYZ))
);

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
       *mag(txz.internalField()/sigma0.internalField()) + SMALL, (1.0/beta.value())
    );

numCyc.correctBoundaryConditions();



/*volScalarField f
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
);*/

f = sigma0/(numCyc * T);

dimensionedScalar cv("cv", (mu * k / gamma) * ((2-2*nu)/ ((1-2*nu)+(2-2*nu)*(n*mu/Kprime))));

Info << "alpha :" << alpha << endl;
Info << "beta :" << beta << endl;
Info << "relDen :" << relDensity << endl;
Info << "cv :" << cv << " m^2/s"<<endl;
Info << " gammaD :" << gammaD<< " N/m^3"<< endl;
Info << " period :" << T<< "s"<<endl;
Info << " G :" << mu<< " N/m^2" <<endl;

runTime.write();


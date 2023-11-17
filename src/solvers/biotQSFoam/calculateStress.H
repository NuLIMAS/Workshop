if (runTime.outputTime())
{

Info<< "Calculating stress field sigma\n" << endl;
volSymmTensorField sigma
(
    IOobject
    (
        "sigma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    sigmaD - I*p
);



Info<< "Calculating  field cv\n" << endl;
volScalarField cv
(
    IOobject
    (
        "cv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    (mu * k / gamma) * ((2-2*nu)/ ((1-2*nu)+(2-2*nu)*(n*mu/Kprime)))
);

tauXZ = sigmaD.component(symmTensor::XZ);
Info<< "Max cv = " << max(cv).value()
    << endl;
runTime.write();
}

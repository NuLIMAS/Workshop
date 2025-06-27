if (runTime.outputTime())
{

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

}

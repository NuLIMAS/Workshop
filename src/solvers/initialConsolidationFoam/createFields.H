Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Calculating stress field sigmaD\n" << endl;
volSymmTensorField sigmaD
(
    IOobject
    (
        "sigmaD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mu*twoSymm(fvc::grad(U)) + lambda*(I*tr(fvc::grad(U)))
);

Info<< "Calculating explicit part of div(sigma) divSigmaExp\n" << endl;
volVectorField divSigmaExp
(
    IOobject
    (
        "divSigmaExp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(sigmaD)
);


volScalarField coordZ
(
    IOobject
    (
        "coordZ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -mesh.C().component(vector::Z),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField dist
(
    IOobject
    (
        "dist",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -mesh.C().component(vector::Z),
    zeroGradientFvPatchScalarField::typeName
);



volScalarField one
(
    IOobject
    (
        "one",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 1),
    zeroGradientFvPatchScalarField::typeName
);



volSymmTensorField epsilon
(
    IOobject
    (
        "epsilon",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    symm(fvc::grad(U))
    
);

volScalarField eps
(
    IOobject
    (
        "eps",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0),
    zeroGradientFvPatchScalarField::typeName
    
);

volScalarField zStrain
(
    IOobject
    (
        "zStrain",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0),
    zeroGradientFvPatchScalarField::typeName
);

/*volScalarField fd
(
    IOobject
    (
        "fd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimLength, 0.15),
    zeroGradientFvPatchScalarField::typeName
);*/

dimensionedScalar a(materialProperties.lookup("a"));
dimensionedScalar b(materialProperties.lookup("b"));
dimensionedScalar c(materialProperties.lookup("c"));

volScalarField q
(
    IOobject
    (
        "q",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    a*log(1+(b*mag(zStrain)))+(c*mag(zStrain))
);



/*volScalarField sigmaA
(
    IOobject
    (
        "sigmaA",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimForce/dimArea, 0.0),
    zeroGradientFvPatchScalarField::typeName
);*/

volScalarField K0
(
    IOobject
    (
        "K0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField gammaDash
(
    IOobject
    (
        "gammaDash",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    (rhos-rhof)*g*(1-n), 
    zeroGradientFvPatchScalarField::typeName
);

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
    gammaD*coordZ*(1+2*k0)/3
);

volScalarField sigmaInit
(
    IOobject
    (
        "sigmaInit",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    ((gammaD* -mesh.C().component(vector::Z)))* 
    (1 +2*K0)/3 
);


volScalarField sigmaA
(
    IOobject
    (
        "sigmaA",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimForce/dimArea, 0.0),
    zeroGradientFvPatchScalarField::typeName

);

if (structure)

{
volScalarField fd
(
    IOobject
    (
        "fd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
   dimensionedScalar(materialProperties.lookup("foundationWidth")),
    zeroGradientFvPatchScalarField::typeName
);
dimensionedScalar QbyBL(materialProperties.lookup("QbyBL"));
dimensionedScalar Q(materialProperties.lookup("Q"));
label patchID = mesh.boundaryMesh().findPatchID("seabedFoundation");
const vectorField& faceCentres = mesh.Cf().boundaryField()[patchID];
forAll(mesh.C(),cellI)
{
    scalar mink= 10.0;
    forAll(faceCentres, facei)
    {
        if ( mag(mesh.C()[cellI] - faceCentres[facei]) < mink)
        {
            mink =mag(mesh.C()[cellI] - faceCentres[facei]) ;
        }
    }
    dist[cellI] = mink;
}

sigmaA.internalField() =QbyBL.value()* (1 - 
    ((pow(dist.internalField(),3.0))/
    (pow(pow(fd.internalField(),2)+ pow(dist.internalField(),2),1.5))));

}


sigmaInit = sigmaInit+ (sigmaA*(1 + 2* K0 )/3 );
sigma0 = sigma0 + (sigmaA*(1 + 2* k0 )/3 );

//sv.internalField() = QbyBL.value()* (1 - 
//    ((pow(dist.internalField(),3.0))/
//    (pow(pow(fd.internalField(),2)+ pow(dist.internalField(),2),1.5))));


if (compactNormalStress)
{
    divSigmaExp -= fvc::laplacian(2*mu + lambda, U, "laplacian(DD,D)");
}
else
{
    divSigmaExp -= fvc::div((2*mu + lambda)*fvc::grad(U), "div(sigmaD)");
}

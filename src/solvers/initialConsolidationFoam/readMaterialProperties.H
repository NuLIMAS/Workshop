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

volScalarField n
(
    IOobject
    (
        "n",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(materialProperties.lookup("n")),
    zeroGradientFvPatchScalarField::typeName
);



vector gravity(0, 0, -9.81);
volVectorField g
(
    IOobject
    (
        "g",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector(
        "g",
        dimensionSet(0, 1, -2 , 0, 0),
        gravity
    ),
    zeroGradientFvPatchTensorField::typeName
 );




dimensionedScalar k0(materialProperties.lookup("k0"));
Switch planeStress(materialProperties.lookup("planeStress"));
dimensionedScalar gamma(materialProperties.lookup("gamma"));

dimensionedScalar Kprime(materialProperties.lookup("K"));


volScalarField rhos
(
    IOobject
    (
        "rhos",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(materialProperties.lookup("rhos")),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField rhof
(
    IOobject
    (
        "rhof",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(materialProperties.lookup("rhof")),
    zeroGradientFvPatchScalarField::typeName
);



volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    n*rhof + (1-n)*rhos
);



volScalarField gammaD
(
    IOobject
    (
        "gammaD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (1-n)*(rhos-rhof)*g.component(vector::Z),
    zeroGradientFvPatchScalarField::typeName
);


    Info<< "Calculating sigma0\n" << endl;
    volScalarField sigma0_static
    (
        IOobject
        (
            "sigma0_static",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (gammaD* mesh.C().component(vector::Z))*
        (1 + 2*k0)/3
    );
    volScalarField one
(
    IOobject
    (
        "one",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 1.0),
    zeroGradientFvPatchScalarField::typeName
);


Switch stones(materialProperties.lookup("stones"));
if(stones)
{
    dimensionedScalar Es(materialProperties.lookup("Es"));
    dimensionedScalar ns(materialProperties.lookup("ns"));
    dimensionedScalar nus(materialProperties.lookup("nus"));
    dimensionedScalar hs(materialProperties.lookup("hs"));
    dimensionedScalar gammaDs(materialProperties.lookup("gammaDs"));


    const cellZone& cz1 = mesh.cellZones()["stone"];

    forAll(cz1, czi)
    {
        E[cz1[czi]] = Es.value();
        nu[cz1[czi]] = nus.value();
        n[cz1[czi]] = ns.value();
        gammaD[cz1[czi]] = gammaDs.value();

    }

forAll (mesh.C(),cellI)
{

if (mesh.C()[cellI].component(vector::Z) > -hs.value() )
{

    sigma0_static[cellI] = gammaDs.value()* one[cellI] * (-mesh.C()[cellI].component(vector::Z)) * ((1.0 + (2.0 * k0.value())) / 3.0);

}
else
{

sigma0_static[cellI] = gammaDs.value() * (hs.value()) * ((1.0 + (2.0 * k0.value())) / 3.0) * one[cellI]
                - gammaD[cellI] * (-hs.value() - mesh.C()[cellI].component(vector::Z)) * ((1.0 + (2.0 * k0.value())) / 3.0);

Info <<gammaD[cellI] * (-hs.value() - mesh.C()[cellI].component(vector::Z)) * ((1.0 + (2.0 * k0.value())) / 3.0) << endl;
            }
}



}


Info << "gammaD" << gammaD << endl;




    volScalarField Erez
    (
        IOobject
        (
            "Erez",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimForce/dimArea, 0.23e6),
        zeroGradientFvPatchScalarField::typeName

    );


    // volScalarField E0
    // (
    //     IOobject
    //     (
    //         "E0",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh,
    //     dimensionedScalar("zero", dimForce/dimArea, 0),
    //     zeroGradientFvPatchScalarField::typeName

    // );

    // E0.internalField() =  6.8* pow(sigma0.internalField(),0.536) + Erez.internalField();
    // E0.correctBoundaryConditions();

    //  E = E0;


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

    // if (runTime.outputTime())
    // {
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

        sigma0 = sigma0_static + mag(tr(sigmaD)) / 3.0;

        // sigmaA = initStress - sigma0;

    // E.internalField() =  E.internalField() * 1/(1 +2.696e3* mag(volStrain.internalField()));

        runTime.writeNow(); // Forces output regardless of controlDict writeInterval rules
        runTime.write();
    // }

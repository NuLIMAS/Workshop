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
            sigmaD
        );

        volSymmTensorField sigmaDev
        (
            IOobject
            (
                "sigmaDev",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            dev(sigma)
        );

        volScalarField sigmaEq
        (
            IOobject
            (
                "sigmaEq",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sqrt((3.0/2.0)*magSqr(dev(sigma)))
        );

        Info<< "Max sigmaEq = " << max(sigmaEq).value()
            << endl;
            
        epsilon = symm(fvc::grad(U));
        eps = -tr(epsilon);
        zStrain = -epsilon.component(symmTensor::ZZ);
        q= a*log(1+(b*eps))+(c*eps);
        
        E == a*b/(1+(b*eps)) + c;
        mu == E/(2.0*(1.0 + nu));
        lambda == nu*E/((1.0 + nu)*(1.0 - 2.0*nu));


        K0 ==(sigmaD.component(symmTensor::XX)/sigmaD.component(symmTensor::ZZ));
        gammaDash =(rhos-rhof)*g*(1-n);
        sigmaInit == ((gammaDash* coordZ))* (1 +2*K0 )/3 ;
        sigma0 = gammaD*coordZ*(1+2*k0)/3;
        sigmaInit = sigmaInit+ (sigmaA*(1 + 2* K0 )/3 );
        sigma0 = sigma0 + (sigmaA*(1 + 2* k0 )/3 );



        
        runTime.write(); 
    }
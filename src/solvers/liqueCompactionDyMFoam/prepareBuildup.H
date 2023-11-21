// Calculate empirical coefficients alpha and beta
volScalarField rd_a = 0.34 * relDensity + 0.084;
volScalarField rd_b = 0.37 * relDensity - 0.46;

// Calculate shear stress components
volScalarField tauXY("tauXY", sigmaD.component(symmTensor::XY));
tauXZ = sigmaD.component(symmTensor::XZ);
volScalarField tauYZ("tauYZ", sigmaD.component(symmTensor::YZ));

Info<< "runTime :" << runTime.time().value() <<endl;

volScalarField txz = SMALL*sigmaD.component(symmTensor::XZ);
if ( runTime.time().value() > 1.0*T.value())
{
     txz = sqrt(mag(mesh.lookupObject<volScalarField>("tauXZPrime2Mean"))*2.0);
}

// Calculate the number of cycles based on shear stress and relative density
numCyc.internalField() =
    Foam::pow
    (
        (1.0/rd_a.internalField())
       *mag(txz.internalField()/sigma0.internalField())
       + SMALL, (1.0/rd_b.internalField())
    );

numCyc.correctBoundaryConditions();

// Calculate source term for pressure accumulation
f= (1-Cl)*sigma0/(numCyc * T);

f.correctBoundaryConditions();

// Solve the governing equation for accumulated pore pressure
fvScalarMatrix pEEqn
(
    fvm::ddt(pE)
    ==
    fvm::laplacian(cv, pE, "laplacian(cv,pE)")
    +  (1.0 - Cl)*f
);
pEEqn.solve();

runTime.write();




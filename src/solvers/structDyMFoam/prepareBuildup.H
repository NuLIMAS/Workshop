
volScalarField rd_a = 0.34 * relDensity + 0.084;
volScalarField rd_b = 0.37 * relDensity - 0.46;

// // Calculate shear stress components
tauXY = sigmaD.component(symmTensor::XY);
tauXZ = sigmaD.component(symmTensor::XZ);
tauYZ = sigmaD.component(symmTensor::YZ);

Info<< "runTime :" << runTime.time().value() <<endl;

volScalarField txz = SMALL*sigmaD.component(symmTensor::XZ);
volScalarField txy = SMALL*sigmaD.component(symmTensor::XY);
volScalarField tyz = SMALL*sigmaD.component(symmTensor::YZ);


if ( runTime.time().value() > 1.0*T.value())
{

    txz = sqrt(mesh.lookupObject<volScalarField>("tauXZPrime2Mean")*2.0);
    txy = sqrt(mesh.lookupObject<volScalarField>("tauXYPrime2Mean")*2.0);
    tyz = sqrt(mesh.lookupObject<volScalarField>("tauYZPrime2Mean")*2.0);
}

// Calculate the number of cycles based on shear stress and relative density
numCycXY.internalField() =
    Foam::pow
    (
        (1.0/rd_a.internalField())
       *mag(txy.internalField()/sigma0.internalField())
       + SMALL, (1.0/rd_b.internalField())
    );
numCycXZ.internalField() =
    Foam::pow
    (
        (1.0/rd_a.internalField())
       *mag(txz.internalField()/sigma0.internalField())
       + SMALL, (1.0/rd_b.internalField())
    );
numCycYZ.internalField() =
    Foam::pow
    (
        (1.0/rd_a.internalField())
       *mag(tyz.internalField()/sigma0.internalField())
       + SMALL, (1.0/rd_b.internalField())
    );
numCycXY.correctBoundaryConditions();
numCycXZ.correctBoundaryConditions();
numCycYZ.correctBoundaryConditions();

// Calculate source term for pressure accumulation
f= (sigma0/(numCycXZ * T)  + sigma0/(numCycYZ* T) +sigma0/(numCycXY * T) );

f.correctBoundaryConditions();

fvScalarMatrix pEEqn
(
    fvm::ddt(pE)
    ==
    fvm::laplacian(cv, pE, "laplacian(cv,pE)")
    +  Cf*f
);
pEEqn.solve();

runTime.write();

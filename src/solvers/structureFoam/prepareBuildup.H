
volScalarField rd_a = 0.34 * relDensity + 0.084;
volScalarField rd_b = 0.37 * relDensity - 0.46;

volScalarField tauXY("tauXY", sigmaD.component(symmTensor::XY));
tauXZ = sigmaD.component(symmTensor::XZ);
volScalarField tauYZ("tauYZ", sigmaD.component(symmTensor::YZ));

Info<< "runTime :" << runTime.time().value() <<endl;

volScalarField txz = SMALL*sigmaD.component(symmTensor::XZ);

if ( runTime.time().value() > T.value())
{
     txz = sqrt(mag(mesh.lookupObject<volScalarField>("tauXZPrime2Mean"))*2.0);
}

numCyc.internalField() =
    Foam::pow
    (
        (1.0/rd_a.internalField())
       *mag(txz.internalField()/sigma0.internalField())
       + SMALL, (1.0/rd_b.internalField())
    );

numCyc.correctBoundaryConditions();

f= sigma0/(numCyc * T);

f.correctBoundaryConditions();

fvScalarMatrix pEEqn
(
    fvm::ddt(pE)
    ==
    fvm::laplacian(cv, pE, "laplacian(cv,pE)")
    +  Cf* f
);
pEEqn.solve();

liqueCriteria.internalField() = pE.internalField()/sigma0.internalField();
liqueCriteria.correctBoundaryConditions();

runTime.write();


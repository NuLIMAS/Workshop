// Calculate shear stress components
volScalarField tauXY("tauXY", sigmaD.component(symmTensor::XY));
tauXZ = sigmaD.component(symmTensor::XZ);
volScalarField tauYZ("tauYZ", sigmaD.component(symmTensor::YZ));

volScalarField txz = SMALL*sigmaD.component(symmTensor::XZ);
if ( runTime.time().value() > T.value())
{
     txz = sqrt(mesh.lookupObject<volScalarField>("tauXZPrime2Mean")*2.0);
}

// Calculate the liqeuafction criteria
liqueCriteria.internalField() = pE.internalField()/sigma0.internalField();
liqueCriteria.correctBoundaryConditions();

// Calculate the number of cycles based on shear stress and relative density
numCyc.internalField() =
    Foam::pow
    (
        (1.0/alpha.value())
       *mag(txz.internalField()/sigma0.internalField()) + SMALL, (1.0/beta.value())
    );

numCyc.correctBoundaryConditions();

// Calculate source term for pressure accumulation
f = sigma0/(numCyc * T);

// Calculate coefficient of consolidation
dimensionedScalar cv("cv", (mu * k / gamma) * ((2-2*nu)/ ((1-2*nu)+(2-2*nu)*(n*mu/Kprime))));

Info << "alpha :" << alpha.value() << endl;
Info << "beta :" << beta.value() << endl;
Info << "cv :" << cv << " m^2/s"<<endl;
runTime.write();


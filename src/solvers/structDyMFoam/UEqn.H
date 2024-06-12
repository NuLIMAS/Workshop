lambda.correctBoundaryConditions();
mu.correctBoundaryConditions();
fvVectorMatrix UEqn
(
    fvm::laplacian(2*mu + lambda, U, "laplacian(DD,U)")
    + divSigmaExp
    - fvm ::d2dt2(rhom,U)
    + Us
    //- rhof*fvc::ddt(Vm)
    //- (rhof/rhom)*fvc::div(phi,Vm)
    ==
    fvc::grad(p)
    - (1-Cl)*fvc::grad(pE)
   - rhodg
);

UEqn.relax();
UResidual = UEqn.solve().initialResidual();
U.relax();

gradU = fvc::grad(U);
sigmaD = mu*twoSymm(gradU) + (lambda*I)*tr(gradU);
divSigmaExp = fvc::div(sigmaD - (2*mu + lambda)*gradU,"div(sigmaD)");
strain = 0.5*twoSymm(fvc::grad(U));
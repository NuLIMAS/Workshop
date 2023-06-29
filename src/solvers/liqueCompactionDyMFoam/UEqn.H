lambda.correctBoundaryConditions();
mu.correctBoundaryConditions();
fvVectorMatrix UEqn
(
    fvm::laplacian(2*mu + lambda, U, "laplacian(DD,U)")
    + divSigmaExp
    ==
    fvc::grad(p)
);

UEqn.relax();
UResidual = UEqn.solve().initialResidual();
U.relax();

gradU = fvc::grad(U);
sigmaD = mu*twoSymm(gradU) + (lambda*I)*tr(gradU);
divSigmaExp = fvc::div(sigmaD - (2*mu + lambda)*gradU,"div(sigmaD)");
strain = 0.5*twoSymm(fvc::grad(U));

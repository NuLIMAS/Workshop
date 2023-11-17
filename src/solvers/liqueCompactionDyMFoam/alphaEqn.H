// Correct boundary conditions for the volume fraction alpha
alpha.correctBoundaryConditions();

// Interpolate alpha to cell faces
alphaF=fvc::interpolate(alpha);

phiAlpha = fvc::interpolate(Cl)*((mesh.Sf() & fvc::interpolate(V))
            + mesh.magSf()*Vdj*(rhoc/fvc::interpolate(rho)));

// Define the alpha transport equation
fvScalarMatrix AlphaEqn
(
    fvm::ddt(alpha)
    + Cl* fvm::div(phiAlpha, alpha)
);

// Solve the alpha transport equation
solve (AlphaEqn);


alpha.correctBoundaryConditions();
alphaF=fvc::interpolate(alpha);

phiAlpha = fvc::interpolate(Cl)*((mesh.Sf() & fvc::interpolate(Vm))
            + mesh.magSf()*Vdj*(rhoc/fvc::interpolate(rho)));

fvScalarMatrix AlphaEqn
(
    Cl*fvm::ddt(alpha)
    +(1-Cl)*fvm::ddt(alpha)
    + Cl* fvm::div(phiAlpha, alpha)
);

solve (AlphaEqn);

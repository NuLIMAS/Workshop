surfaceScalarField phiVdj ("phiVdj",fvc::interpolate(alpha/(1-alpha)*rhoc*rhos/rho));

fvVectorMatrix VEqn
(
      Cl * fvm::ddt(rhom, V)
    + Cl * fvm::div(phi, V)
    - Cl * fvm::laplacian(visf, V)
    + Cl * fvc::div
    (
        phiVdj* (mesh.Sf()*Vdj*Vdj)
    )
    + (1.0 - Cl)*(rho)*fvm::Sp((mag(g)/k), V)
    - (1.0 - Cl)*(1-Cs)* fvc:: grad(pE)
);
VEqn.relax();

// Momentum solution
if (piso.momentumPredictor())
{
    solve(VEqn == - fvc::grad(p));
}
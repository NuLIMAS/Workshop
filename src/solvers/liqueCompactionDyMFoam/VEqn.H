surfaceScalarField phiVdj ("phiVdj",fvc::interpolate(alpha/(1-alpha)*rhoc*rhos/rho));

fvVectorMatrix VmEqn
(
      Cl * fvm::ddt(rhom, Vm)
    + Cl * fvm::div(phi, Vm)
    - Cl * fvm::laplacian(visf, Vm)
    + Cl * fvc::div
    (
        phiVdj* (mesh.Sf()*Vdj*Vdj)
    )
    + (1.0 - Cl)*(rho)*fvm::Sp((mag(g)/k), Vm)
    - (1.0 - Cl)*(1-Cs)* fvc:: grad(pE)
);
VmEqn.relax();

// Momentum solution
if (piso.momentumPredictor())
{
    solve(VmEqn == - fvc::grad(p));
}
volScalarField rAVm("rAVm",1.0/VmEqn.A());

surfaceScalarField rAVmf_temp
(
    "(rho*(1|A(Vmf)))",
    fvc::interpolate(rAVm, "myInterpolation")
);

surfaceScalarField rhoF
(
    "rhoF",
    fvc::interpolate(rho,  "myInterpolation")
);
surfaceScalarField rAVmf ( "rAVmf", rAVmf_temp*rhoF);

Vm = rAVm*VmEqn.H();

phi = rhoF*(fvc::interpolate(Vm,"myInterpolation") & mesh.Sf());
surfaceVectorField  rhosc ("rhosc",(rhoF - rhoc)*g);
phi += rhosc*rAVmf& mesh.Sf();

// Non-orthogonal pressure corrector loop
while (piso.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        - (1-st)*(1.0 - Cl)*(rho*n/Kprime)*fvm::ddt(p)
        + (1-st)*(1.0 - Cl)*(1-Cs)*(rho*n/Kprime)*fvc::ddt(pE)
        - (1-st)*(1.0 - Cl)*(rho)*fvc::div(fvc::ddt(U))
        + fvm::laplacian(rAVmf, p ,"laplacian(rAVmf,p)")
        - Cl * fvc::ddt(rhom)
        ==
        fvc::div(phi)
    );
    pEqn.relax();

    pResidual=pEqn.solve().initialResidual();

    if (piso.finalNonOrthogonalIter())
    {
        phi -= pEqn.flux();
    }
}

// Calculate the continuity error
contErr = fvc::div(phi)
    + (1-st)*(1.0 - Cl)*(rho*n/Kprime)*fvc::ddt(p)
    - (1-st)*(1.0 - Cl)*(1-Cs)*(rho*n/Kprime)*fvc::ddt(pE)
    + (1-st)*(1.0 - Cl)*(rho)*fvc::div(fvc::ddt(U))
    + Cl * fvc::ddt(rhom);

sumLocalContErr = runTime.deltaT().value()*
    mag(contErr)().weightedAverage(mesh.V()).value();

globalContErr = runTime.deltaT().value()*
    contErr.weightedAverage(mesh.V()).value();

cumulativeContErr += globalContErr;

Info<< "time step continuity errors : sum local = " << sumLocalContErr
    << ", global = " << globalContErr
    << ", cumulative = " << cumulativeContErr
    << endl;

p.relax();
surfaceScalarField phiByRho = (phi)/rhoF;
Vm = fvc::reconstruct(phiByRho);
Vm.correctBoundaryConditions();

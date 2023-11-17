// Define and interpolate rAV for non-orthogonal pressure correction
volScalarField rAV("rAV",1.0/VEqn.A());

surfaceScalarField rAVf_temp
(
    "(rho*(1|A(Vf)))",
    fvc::interpolate(rAV, "myInterpolation")
);

surfaceScalarField rhoF
(
    "rhoF",
    fvc::interpolate(rho,  "myInterpolation")
);
surfaceScalarField rAVf ( "rAVf", rAVf_temp*rhoF);

V = rAV*VEqn.H();

// Calculate fluxes for non-orthogonal pressure correction
phi = rhoF*(fvc::interpolate(V,"myInterpolation") & mesh.Sf());
surfaceVectorField  rhosc ("rhosc",(rhoF - rhoc)*g);
phi += rhosc*rAVf& mesh.Sf();

// Non-orthogonal pressure corrector loop
while (piso.correctNonOrthogonal())
{
    // Define pressure equation
    fvScalarMatrix pEqn
    (
        - (1-st)*(1.0 - Cl)*(rho*n/Kprime)*fvm::ddt(p)
        + (1-st)*(1.0 - Cl)*(1-Cs)*(rho*n/Kprime)*fvc::ddt(pE)
        - (1-st)*(1.0 - Cl)*(rho)*fvc::div(fvc::ddt(U))
        + fvm::laplacian(rAVf, p ,"laplacian(rAVmf,p)")
        - Cl * fvc::ddt(rhom)
        ==
        fvc::div(phi)
    );
    pEqn.relax();

    // Solve the pressure equation
    pResidual=pEqn.solve().initialResidual();

    // Update phi for continuity errors
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

// Relax pressure field and update velocity field
p.relax();
surfaceScalarField phiByRho = (phi)/rhoF;
V = fvc::reconstruct(phiByRho);
V.correctBoundaryConditions();

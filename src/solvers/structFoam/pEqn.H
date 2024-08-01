
fvScalarMatrix pEqn
(
    (1-st)*(1/Dp2)*fvm::ddt(p)
    ==
     fvm::laplacian(Dp3, p)
    - (1-st)*  fvc::div(fvc::ddt(U))

);
pEqn.relax();
pResidual = pEqn.solve().initialResidual();
p.relax();

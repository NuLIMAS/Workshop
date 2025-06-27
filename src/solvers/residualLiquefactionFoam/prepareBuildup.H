
volScalarField tauXY("tauXY", sigmaD.component(symmTensor::XY));
tauXZ = sigmaD.component(symmTensor::XZ);
volScalarField tauYZ("tauYZ", sigmaD.component(symmTensor::YZ));

Info<< "runTime :" << runTime.time().value() <<endl;

volScalarField txz = SMALL*sigmaD.component(symmTensor::XZ);

volScalarField rd_a = 0.34 * relDensity + 0.084;
volScalarField rd_b = 0.37 * relDensity - 0.46;

if ( runTime.time().value() > 1.0*T.value())
{
    txz = sqrt(mesh.lookupObject<volScalarField>("tauXZPrime2Mean")*2.0);
}

liqueCriteria.internalField() = pE.internalField()/sigma0.internalField();
liqueCriteria.correctBoundaryConditions();



numCyc.internalField() =
Foam::pow
(
    (1.0/rd_a.internalField())
    *mag(txz.internalField()/sigma0.internalField()) + SMALL, (1.0/rd_b.internalField())
);

numCyc.correctBoundaryConditions();

f = sigma0/(numCyc * T);

scalarField& ClI = Cl.internalField();

forAll (ClI,cellI)
{
    if (  liqueCriteria[cellI] >= 1.0 )
    {
        ClI[cellI] =1;
        liqueFlag[cellI]=0;
    }
    else
    {
        ClI[cellI] =0;
        liqueFlag[cellI]=1;
    }
}

Cl.correctBoundaryConditions();
liqueFlag.correctBoundaryConditions();

volScalarField cv("cv", (mu * k / gamma) * ((2-2*nu)/ ((1-2*nu)+(2-2*nu)*(n*mu/Kprime))));

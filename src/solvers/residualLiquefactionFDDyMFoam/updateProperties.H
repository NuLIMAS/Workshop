

scalarField& Cli = Cl.internalField();

forAll (Cli,cellI)
{

if (  Cl[cellI] >= 1.0 )
{
    //Info << " Liquefaction occurs" << endl;

    //E[cellI]=2e7;
    //k[cellI]=0.004;
    //n[cellI]=0.7;

}
}

//Dp4=k/g;
//Dp =k/gamma*Kprime/n;
//Dp2=Kprime/n;
//Dp3=k/gamma;

mu = E/(2.0*(1.0 + nu));
lambda = nu*E/((1.0 + nu)*(1.0 - 2.0*nu));




Switch planeStress(materialProperties.lookup("planeStress"));

if (planeStress)
{
    Info<< "Plane Stress\n" << endl;
    lambda = nu*E/((1.0 + nu)*(1.0 - nu));
}
else
{
    Info<< "Plane Strain\n" << endl;
}
    

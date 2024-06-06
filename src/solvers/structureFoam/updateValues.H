/*forAll (mesh.C(),cellI)
{
    if ( mesh.C().internalField()[cellI].component(vector::Z) > 0 )
    {
        //rhom[cellI] =2656;
        E[cellI] = 1e10;
        Cp[cellI] = 0;
        rhom[cellI] =1395;
        Cf[cellI] =0;
        //E[cellI] = 1e9;
        //Us[cellI] = vector(0,0,-25000);
    }



}*/

/*
forAll (mesh.C(),cellI)
{
    if ( mesh.C().internalField()[cellI].component(vector::Z) > 0 )
    {
        //rhom[cellI] =2656;
        E[cellI] = 1e9;
        Cp[cellI] = 0;
        rhom[cellI] =2656;
        nu[cellI] = 0.15;
        k[cellI] = 1e-7;
        n[cellI] = 0.25;
        //E[cellI] = 1e9;
        //Us[cellI] = vector(0,0,-2500);

    }
}*/


Switch stationary(materialProperties.lookup("stationary"));
label boundaryID3 = mesh.boundaryMesh().findPatchID("foundationSeabed");
//dimensionedScalar us(materialProperties.lookup("Us"));

int tp = ceil( runTime.time().value() / T.value()) ;

const cellZone& cz1 = mesh.cellZones()["structureZones"];
forAll(cz1, czi)
{
    E[cz1[czi]] = 1e10;
    rhom[cz1[czi]] =rhoa.value();
    rho[cz1[czi]] =rhoa.value();
    Cf[cz1[czi]] =0;
    Cl[cz1[czi]] =0;
    //E[cellI] = 1e10;
    //rhom[cellI] =1450;
    //Cf[cellI] =0;
    //Cl[cellI] =0;
}
rho.boundaryField()[boundaryID3] = rhoa.value();
rhom.boundaryField()[boundaryID3] = rhoa.value();
const cellZone& cz2 = mesh.cellZones()["ropeZones1"];
forAll(cz2, czi)
{
    if (tp %2 == 0)
    {
        Us[cz2[czi]] = Ur1.value();
    }
    else
    {
        Us[cz2[czi]] = Ur2.value();
    }
}
const cellZone& cz3 = mesh.cellZones()["ropeZones2"];
forAll(cz3, czi)
{
    if (tp %2 == 0)
    {
        Us[cz3[czi]] = Ur2.value();
    }
    else
    {
        Us[cz3[czi]] = Ur1.value();
    }
}

/*if (stationary)
{
    forAll (mesh.C(),cellI)
    {
        if ( mesh.C().internalField()[cellI].component(vector::Z) > 0 )
        {
            //rhom[cellI] =2656;
            E[cellI] = 1e10;
            //Cp[cellI] = 0;
            rhom[cellI] =1395;
            Cf[cellI] =0;
            //E[cellI] = 1e9;
            //Us[cellI] = vector(0,0,-25000);
        }
    }
}
else
{
    //Info << "Moving foundation" << endl;
    forAll (mesh.C(),cellI)
    {
        if  (mesh.C().internalField()[cellI].component(vector::Z) >  0.035 &&   // 0.028 && //0.009
            mesh.C().internalField()[cellI].component(vector::X) > 1.02  &&//0.35 &&
            mesh.C().internalField()[cellI].component(vector::X) <  1.03 //0.36     //0.38  0.353
            )
        {
            if (tp %2 == 0)
            {
                Us[cellI] = vector(0,0,us.value());
            }
            else
            {
                Us[cellI] = vector(0,0,0);
            }
        }
        if  (mesh.C().internalField()[cellI].component(vector::Z) > 0.035 &&
            mesh.C().internalField()[cellI].component(vector::X) > 1.29 &&//0.64 &&  //0.62  && 0.647 
             mesh.C().internalField()[cellI].component(vector::X) < 1.3 //0.65     
            )
        {
            if (tp %2 == 0)
            {
                Us[cellI] = vector(0,0,0);
            }
            else
            {
                Us[cellI] = vector(0,0,us.value());
            }
        }
        //if ( mesh.C().internalField()[cellI].component(vector::Z) > 0.035 )
        //{
        //    Us[cellI] = vector(0,0,125000);
        //}

        if ( mesh.C().internalField()[cellI].component(vector::Z) > 0 )
        {
            E[cellI] = 1e10;
            rhom[cellI] =rhoa.value();
            Cf[cellI] =0;
        }

    }
}*/


E.correctBoundaryConditions();
rhom.correctBoundaryConditions();
Cf.correctBoundaryConditions();


// Update Shear modulus and Lames constant
mu = E/(2.0*(1.0 + nu));
lambda = nu*E/((1.0 + nu)*(1.0 - 2.0*nu));


if (planeStress)
{
    lambda = nu*E/((1.0 + nu)*(1.0 - nu));
}

// Update Coefficient of consolidation
cv = (mu * k / gamma) * ((2-2*nu)/ ((1-2*nu)+(2-2*nu)*(n*mu/Kprime)));
cv.correctBoundaryConditions();
//volScalarField Dp  =  (k/gamma*Kprime/n);
//volScalarField Dp2 =  (Kprime/n);
//volScalarField Dp2 ("Dp2", Kprime/n);
//dimensionedScalar

volScalarField Dp2
(
    "Dp2",
    (Kprime/n)
);

volScalarField Dp3
(
    "Dp3",
    (k/gamma)
);
Dp2.correctBoundaryConditions();
Dp3.correctBoundaryConditions();
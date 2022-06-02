

Info<< "Reading wave properties\n" << endl;

IOdictionary waveProperties
(
    IOobject
    (
        "waveProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);




dimensionedScalar e = n/(1.0 -n);
dimensionedScalar k0(materialProperties.lookup("k0"));
dimensionedScalar relDensity(materialProperties.lookup("relDensity"));
//dimensionedScalar alpha(materialProperties.lookup("alpha"));
//dimensionedScalar beta(materialProperties.lookup("beta"));
dimensionedScalar alpha = 0.34 * relDensity + 0.084;
dimensionedScalar beta = 0.37 * relDensity - 0.46;

dimensionedScalar T("T",dimensionSet(0, 0, 1, 0, 0), waveProperties.lookup("period") );

volScalarField tauXY("tauXY", sigmaD.component(symmTensor::XY));
tauXZ = sigmaD.component(symmTensor::XZ);
volScalarField tauYZ("tauYZ", sigmaD.component(symmTensor::YZ));

dimensionedScalar gammaD(materialProperties.lookup("gammaD"));
Info<< "runTime :" << runTime.time().value() <<endl;

volScalarField txz = SMALL*sigmaD.component(symmTensor::XZ);

volScalarField tau_xyz = sqrt(tauXY*tauXY + tauXZ*tauXZ + tauYZ*tauYZ );
if ( runTime.time().value() > T.value())
{
     // txz = sqrt(sigmaD.component(symmTensor::XZ)*sigmaD.component(symmTensor::XZ));
     txz = sqrt(mesh.lookupObject<volScalarField>("tauXZPrime2Mean")*2.0);
}

volScalarField sigma0
(
    IOobject
    (
        "sigma0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    gammaD*(max(mesh.C().component(vector::Z)) -    mesh.C().component(vector::Z))*((1.0+(2.0*k0))/3.0)
);


volScalarField tauEq
(
    IOobject
    (
        "tauEq",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    sqrt((tauXY*tauXY)+(tauXZ*tauXZ)+(tauYZ*tauYZ))
);

liqueCriteria.internalField() = pE.internalField()/sigma0.internalField();
liqueCriteria.correctBoundaryConditions();

volScalarField numCyc
(
    IOobject
    (
        "numCyc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

numCyc.internalField() =
    Foam::pow
    (
        (1.0/alpha.value())
       *mag(txz.internalField()/sigma0.internalField()) + SMALL, (1.0/beta.value())
    );

numCyc.correctBoundaryConditions();

forAll (mesh.C(),cellI)
{
     if ( liqueCriteria[cellI] > 1.0 )
     {
         Cs[cellI] =0;
         liqueFlag[cellI] =0;          
     }
     
     if ( mesh.C().internalField()[cellI].component(vector::Z) > -0.015 &&
            mesh.C().internalField()[cellI].component(vector::X) > 0.3  &&  
            mesh.C().internalField()[cellI].component(vector::X) < 0.6  )
     {
         Cs[cellI] =0; 
     }
     
     if ( mesh.C().internalField()[cellI].component(vector::Z) > -0.03 &&
           mesh.C().internalField()[cellI].component(vector::X) > 0.25 && 
           mesh.C().internalField()[cellI].component(vector::X) < 0.30 &&
           mesh.C().internalField()[cellI].component(vector::X) > 0.6  &&     
           mesh.C().internalField()[cellI].component(vector::X) < 0.65  )
     {
         Cs[cellI] =0; 
     }

}


volScalarField f
(
    IOobject
    (
        "f",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    Cs*sigma0/(numCyc * T)
);

const scalarField& iLiq = liqueFlag.internalField();
const labelList& nei = mesh.neighbour();
const labelList& own = mesh.owner();
DynamicList<label> newFace(nei.size());
forAll(nei, facei)
{
    // Adress volField neighbour side
    iLiq[nei[facei]]; 
    iLiq[own[facei]];
     if  (iLiq[nei[facei]] == 1 &&  iLiq[own[facei]] == 0)
     {
         newFace.append(facei);
     }
     else if  (iLiq[nei[facei]] == 0 &&  iLiq[own[facei]] == 1)
     {
          newFace.append(facei);
     }
            
}


dimensionedScalar cv("cv", (mu * k / gamma) * ((2-2*nu)/ ((1-2*nu)+(2-2*nu)*(n*mu/Kprime))));

//dimensionedScalar cv("cv", (mu * k / gamma) * (1.0/(1-2*nu)));

Info << "alpha :" << alpha << endl;
Info << "beta :" << beta << endl;
Info << "relDen :" << relDensity << endl;
Info << "cv :" << cv << " m^2/s"<<endl;
Info << " gammaD :" << gammaD<< " N/m^3"<< endl;
Info << " period :" << T<< "s"<<endl;
Info << " G :" << mu<< " N/m^2" <<endl;







runTime.write();


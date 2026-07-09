Switch stones(materialProperties.lookupOrDefault<bool>("stones",0));
if (stones)
{
    dimensionedScalar Es(materialProperties.lookup("Es"));
    dimensionedScalar ns(materialProperties.lookup("ns"));
    dimensionedScalar ks(materialProperties.lookup("ks"));
    dimensionedScalar relDensitys(materialProperties.lookup("relDensitys"));


    const cellZone& cz1 = mesh.cellZones()["stoneZone"];
    forAll(cz1, czi)
    {

        E[cz1[czi]] = Es.value();
        n[cz1[czi]] = ns.value();
        alpha[cz1[czi]] = 1.0 - ns.value();
        k[cz1[czi]] = ks.value();
        Cf[cz1[czi]] = 0;
        relDensity[cz1[czi]]= relDensitys.value();
    }
}

surfaceScalarField alphaF
(
    "alphaF",
    fvc::interpolate(alpha)
);

// Bound alphaF to prevent log of non-positive number and other math errors
alphaF = min(1.0 - SMALL, max(0.0, alphaF));
Vdj = V0 & mesh.Sf()/mesh.magSf()*pow(1-alphaF,q);


label boundaryID1 = mesh.boundaryMesh().findPatchID("sol");
label boundaryID2 = mesh.boundaryMesh().findPatchID("liq");
Vdj.boundaryField()[boundaryID1] = 0;
Vdj.boundaryField()[boundaryID2] = 0;
const fvPatchScalarField& patchUsol = Cl.boundaryField()[boundaryID1];
const fvPatchScalarField& patchUliq = Cl.boundaryField()[boundaryID2];
const labelList Clsol=patchUsol.patch().faceCells();
const labelList Clliq=patchUliq.patch().faceCells();

// Update viscosity and density

//visf = muf*(1+2.5*alphaF+10.05*alphaF*alphaF + 0.00273*exp(16.6*alphaF));
visf = muf + C* log(max(1.0 - alphaF, SMALL));
rhom = (1.0 - alpha) * rhoc + alpha * rhos;
rho = ((1.0 - alpha) * rhoc + alpha * rhos) * Cl + (1-Cl) * rhoc;

//Update density, Cl at the boundary
forAll(patchUliq, i)
{
    Cl.boundaryField()[boundaryID2][i] = 1;
    Cl.boundaryField()[boundaryID1][i] = 0;
    rho.boundaryField()[boundaryID2][i] = rhoc.value();
    rho.boundaryField()[boundaryID1][i] = rhoc.value();
}

E.correctBoundaryConditions();
nu.correctBoundaryConditions();
relDensity.correctBoundaryConditions();


Cl.correctBoundaryConditions();
n.correctBoundaryConditions();
//Cs.correctBoundaryConditions();
//rho.correctBoundaryConditions();
Vdj.correctBoundaryConditions();

mu = E/(2.0*(1.0 + nu));
lambda = nu*E/((1.0 + nu)*(1.0 - 2.0*nu));

Switch planeStress(materialProperties.lookup("planeStress"));

if (planeStress)
{
//    Info<< "Plane Stress\n" << endl;
    lambda = nu*E/((1.0 + nu)*(1.0 - nu));
}

// Update Coefficient of consolidation
cv = (mu * k / gamma) * ((2-2*nu)/ ((1-2*nu)+(2-2*nu)*(n*mu/Kprime)));
cv.correctBoundaryConditions();
#if DEBUG
Info << "------ Material properties  --------" << endl;
Info << "Elastic modulus  :   max: " <<  max(E).value() <<  " Pa min: " << min(E).value() << " Pa" << endl;
Info << "Shear modulus    :   max: " <<  max(mu).value() <<  " Pa min: " << min(mu).value() << " Pa" << endl;
Info << "Lames constant   :   max: " <<  max(lambda).value() <<  " Pa min: " << min(lambda).value() << " Pa" << endl;
Info << "Porosity         :   max: " <<  max(n).value() <<  " min: " << min(n).value() << "" << endl;
Info << "Poisson ratio    :   max: " <<  max(nu).value() <<  " min: " << min(nu).value() << "" << endl;
Info << "Permeability     :   max: " <<  max(k).value() <<  " m/s min: " << min(k).value() << " m/s" << endl;
Info << "Density          :   max: " <<  max(rho).value() <<  " kg/m^3 min: " << min(rho).value() << " kg/m^3" << endl;

Info << "--------------------------------------------------" << endl;

Info << "---- Pressurebuildup properties ----" << endl;
Info << "relativeDensity  :   max: " <<  max(relDensity).value() <<  " min: " << min(relDensity).value() << " " << endl;
Info << "Coefficient of lateral Earth pressure  : " <<  k0.value()  << endl;
Info << "Coefficient of consolidation  :   max: " <<  max(cv).value() <<  " m^2/s min: " << min(cv).value() << " m^2/s" << endl;
Info << "---------------------------------------------------" << endl;
#endif



// 1. Force parallel synchronization of cell states across processor boundaries
Cl.correctBoundaryConditions();

// --- 2. MISMATCH-ONLY OVERRIDE FOR ALPHA PROC PATCHES ---
forAll(alpha.boundaryField(), patchi)
{
    if (alpha.mesh().boundary()[patchi].coupled())
    {
        fvPatchScalarField& alphaPatch = alpha.boundaryField()[patchi];

        const labelList& faceCells = alpha.mesh().boundary()[patchi].faceCells();
        scalarField neighborCl = alpha.mesh().boundary()[patchi].lookupPatchField<volScalarField, scalarField>("Cl");
        scalarField neighborAlpha = alpha.mesh().boundary()[patchi].lookupPatchField<volScalarField, scalarField>("alpha");

        forAll(faceCells, facei)
        {
            label localCellID = faceCells[facei];

            // Condition: Local is Solid (0) AND Neighbor is Liquid (1)
            // OR Local is Liquid (1) AND Neighbor is Solid (0)
            if ((Cl[localCellID] == 0 && neighborCl[facei] == 1) || 
                (Cl[localCellID] == 1 && neighborCl[facei] == 0))
            {
                alphaPatch[facei] = 0.0; // Enforce fixedValue 0 at the phase-change interface
            }
            else
            {
                // No mismatch (both solid or both liquid): map across cores normally
                alphaPatch[facei] = neighborAlpha[facei];
            }
        }
    }
}

// --- 3. MISMATCH-ONLY OVERRIDE FOR PE (EXCESS PORE PRESSURE) PROC PATCHES ---
forAll(pE.boundaryField(), patchi)
{
    if (pE.mesh().boundary()[patchi].coupled())
    {
        Info << "Processing coupled patch for pE: " << pE.mesh().boundary()[patchi].name() << endl;
        fvPatchScalarField& pEPatch = pE.boundaryField()[patchi];

        const labelList& faceCells = pE.mesh().boundary()[patchi].faceCells();
        scalarField neighborCl = pE.mesh().boundary()[patchi].lookupPatchField<volScalarField, scalarField>("Cl");
        scalarField neighborPE = pE.mesh().boundary()[patchi].lookupPatchField<volScalarField, scalarField>("pE");

        forAll(faceCells, facei)
        {
            label localCellID = faceCells[facei];

            // Mismatch Scenario A: Local is Liquid (1) AND Neighbor is Solid (0)
            if (Cl[localCellID] == 1 && neighborCl[facei] == 0)
            {
                Info << "Mismatch detected at face " << facei << ": Local is Liquid, Neighbor is Solid." << endl;
                pEPatch[facei] = 0.0; // Acts like the 'liq' side: fixedValue 0
            }
            // Mismatch Scenario B: Local is Solid (0) AND Neighbor is Liquid (1)
            else if (Cl[localCellID] == 0 && neighborCl[facei] == 1)
            {
                Info << "Mismatch detected at face " << facei << ": Local is Solid, Neighbor is Liquid." << endl;
                pEPatch[facei] = pE[localCellID]; // Acts like the 'sol' side: zeroGradient
            }
            else
            {
                // No mismatch (both solid or both liquid): map across cores normally
                pEPatch[facei] = neighborPE[facei];
            }
        }
    }
}

// --- 4. PRINT MISMATCH PE VALUES ON PROC 0 (for debugging) ---
if (Pstream::myProcNo() == 0)
{
    Info << "--- Mismatch pE values on proc 0 ---" << endl;
    forAll(pE.boundaryField(), patchi)
    {
        if (pE.mesh().boundary()[patchi].coupled())
        {
            const fvPatchScalarField& pEPatch = pE.boundaryField()[patchi];
            const labelList& faceCells = pE.mesh().boundary()[patchi].faceCells();
            scalarField neighborCl = pE.mesh().boundary()[patchi].lookupPatchField<volScalarField, scalarField>("Cl");

            forAll(faceCells, facei)
            {
                label localCellID = faceCells[facei];
                if ((Cl[localCellID] == 0 && neighborCl[facei] == 1) ||
                    (Cl[localCellID] == 1 && neighborCl[facei] == 0))
                {
                    Info << "Patch " << pE.mesh().boundary()[patchi].name()
                         << " face " << facei << " pE=" << pEPatch[facei] << endl;
                }
            }
        }
    }
    Info << "------------------------------------" << endl;
}

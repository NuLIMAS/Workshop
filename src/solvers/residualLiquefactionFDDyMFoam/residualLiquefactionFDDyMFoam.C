/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    biotFoam is a modified version of solidDisplacementFoam

Description
    Transient segregated finite-volume solver of the Biot quasi-steady consolidation equations
    for linear-elastic, small-strain deformation of a solid skeleton coupled with
    pore water flow and pressure governed by Darcy's law.

    Simple linear elasticity structural analysis code.
    Solves for the displacement vector field D and the pore water pressure p.
    Also generating the stress tensor field sigma.
    
    An additional laplacian equation is solved in order to calculate the accumulated 
    pore pressure pE.

    Post-liquefaction, the liquefied region are solved using N-S equation (PISO loop). 
Author
    R. Shanmugasundaram, Wikki GmbH
    H. Rusche, Wikki GmbH

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "zeroGradientFvPatchFields.H"

#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "readMaterialProperties.H"
    #include "readPoroElasticControls.H"
    pisoControl piso(mesh);
    #include "createFields.H"
    #include "initContinuityErrs.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating displacement field\n" << endl;
            scalar sumLocalphiV = 0;


    while (runTime.loop())
    {
        #include "updateProperties.H"
        Info<< "Iteration: " << runTime.value() << nl << endl;
        bool meshChanged = mesh.update();
        mesh.checkMesh(true);

        HashTable<volScalarField*> vsFields(mesh.lookupClass<volScalarField>());

        for
        (
            HashTable<volScalarField*>::
            iterator fieldIter = vsFields.begin();
            fieldIter != vsFields.end();
            ++fieldIter
        )
        {
            volScalarField& field = const_cast<volScalarField&>(*fieldIter());
            field.correctBoundaryConditions();
        }

        HashTable<volVectorField*> vvFields(mesh.lookupClass<volVectorField>());

        for
        (
            HashTable<volVectorField*>::
            iterator fieldIter = vvFields.begin();
            fieldIter != vvFields.end();
            ++fieldIter
        )
        {
            volVectorField& field = const_cast<volVectorField&>(*fieldIter());
            field.correctBoundaryConditions();
        }

        HashTable<volSymmTensorField*> vstFields(mesh.lookupClass<volSymmTensorField>());

        for
        (
            HashTable<volSymmTensorField*>::
            iterator fieldIter = vstFields.begin();
            fieldIter != vstFields.end();
            ++fieldIter
        )
        {
            volSymmTensorField& field = const_cast<volSymmTensorField&>(*fieldIter());
            field.correctBoundaryConditions();
        }

        #include "readPoroElasticControls.H"
        


        int iCorr = 0;
        scalar UResidual = 1.0e10;
        scalar pResidual = 1.0e10;
        scalar UfResidual = 1.0e10;
        scalar residual = 1.0e10;
        

        do
        {
            
            U.storePrevIter();
            p.storePrevIter();
            Uf.storePrevIter();
            //Info << Dp4<< endl;
            fvVectorMatrix UfEqn
            (
                Cl*(rhof/n)*fvm::ddt(Uf)
              + Cl*(rhof/n)*fvm::div(phi, Uf)
              + (1.0-Cl)*rhof*fvm::Sp(1.0/Dp4, Uf)
              - Cl*rhof*fvm::laplacian(vis, Uf)
            );
            //UfEqn.relax();

            // Momentum solution
            if (piso.momentumPredictor())
            {
                solve(UfEqn == -fvc::grad(p));
            } 


            //UfEqn.clear();
            while (piso.correct())
            {
                volScalarField rAUf(1.0/UfEqn.A());
                Uf = rAUf*UfEqn.H();
                phi = (fvc::interpolate(Uf) & mesh.Sf());
                // Calculate under-relaxation consistent flux
                adjustPhi(phi, Uf, p);
                // Non-orthogonal pressure corrector loop
                while (piso.correctNonOrthogonal())
                {    
                    fvScalarMatrix pEqn
                    (
                       - (1.0-Cl)*(1/Dp2)*fvm::ddt(p)
                       + fvm::laplacian(rAUf, p ,"laplacian(DD,p)")
                       - (1.0-Cl)*fvc::div(fvc::ddt(U))
                       == 
                       fvc::div(phi) //-Cl*fvm::laplacian(rAUf, p) 
                    );
                    pResidual=pEqn.solve().initialResidual();    


                    if (piso.finalNonOrthogonalIter())
                    {
                        phi -= pEqn.flux();
                    }
               }

               {
               volScalarField contErr = fvc::div(phi)
                   + (1.0-Cl)*(1/Dp2)*fvc::ddt(p)
                   + (1.0-Cl)*fvc::div(fvc::ddt(U));
    
               sumLocalContErr = runTime.deltaT().value()*
                   mag(contErr)().weightedAverage(mesh.V()).value();

               globalContErr = runTime.deltaT().value()*
                   contErr.weightedAverage(mesh.V()).value();

               cumulativeContErr += globalContErr;

               Info<< "time step continuity errors : sum local = " << sumLocalContErr
                   << ", global = " << globalContErr
                   << ", cumulative = " << cumulativeContErr
                   << endl;
               }


               p.relax();
               Uf -= rAUf*fvc::grad(p);
               Uf.correctBoundaryConditions(); 
            }   
                   
            fvVectorMatrix UEqn
            (
               fvm::laplacian(2*mu + lambda, U, "laplacian(DD,U)") 
              + divSigmaExp 
              //-Cl*(rhof)*fvc::ddt(Uf)
              //-Cl*(rhof)*fvc::div(phi, Uf)
              == 
              fvc::grad(p)
            );
            UEqn.relax();
            UResidual = UEqn.solve().initialResidual();
            U.relax();

            gradU = fvc::grad(U);
            sigmaD = mu*twoSymm(gradU) + (lambda*I)*tr(gradU);
            //volScalarField tauXZ=sigmaD.component(symmTensor::XZ);
            divSigmaExp = fvc::div(sigmaD - (2*mu + lambda)*gradU,"div(sigmaD)");

            residual = max(pResidual,UResidual);
  
  
        } while (residual > convergenceTolerance && ++iCorr < nCorr);
        
       
        volScalarField divPhi(
            "divPhi",
            fvc::div(phi)
            + (1.0-Cl)*(1/Dp2)*fvc::ddt(p)
            + (1.0-Cl)*fvc::div(fvc::ddt(U))
        );

        //volScalarField divU = fvc::div(Uf);
        //volScalarField ddtStrain = (1.0-Cl)*fvc::div(fvc::ddt(U));
        //volScalarField ddtp = (1.0-Cl)*(1/Dp2)*fvc::ddt(p);
        //volScalarField divUf(
        //    "divUf",
        //    fvc::div(Uf)
        //);
        //
        //volScalarField ddtStrain(
        //    "ddtStrain",
        //    (1.0-Cl)*fvc::div(fvc::ddt(U))
        //);           
        //  
        //volScalarField ddtp(
        //    "ddtp",
        //    (1.0-Cl)*(1/Dp2)*fvc::ddt(p)
        //);
        // 
        //divPhi.write();
        //divUf.write();
        //ddtStrain.write();
        //ddtp.write();

        Info << "Number of iterations:" << iCorr <<endl;    
        #include "calculateStress.H"
        #include "prepareBuildup.H"

        fvScalarMatrix pEEqn
        (
            fvm::ddt(pE) 
            ==
            fvm::laplacian(cv, pE, "laplacian(cv,pE)")
            + f
        );
        pEEqn.solve();


        
       //std::ofstream file;
       //file.open ("results.txt", std::ofstream::out | std::ofstream::app);
       //file << runTime.timeName() << " " << sumLocalContErr  << std::endl << "\n";
       //file.close();
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

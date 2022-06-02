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
Author
    R. Shanmugasundaram, Wikki GmbH
    H. Rusche, Wikki GmbH
    Johan Roenby, DHI Water & Environment

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "readMaterialProperties.H"
    #include "readPoroElasticControls.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating displacement field\n" << endl;

    while (runTime.loop())
    {
        Info<< "Iteration: " << runTime.value() << nl << endl;

        #include "readPoroElasticControls.H"

        int iCorr = 0;
        scalar UResidual = 1.0e10;
        scalar pResidual = 1.0e10;
        scalar residual = 1.0e10;

        do
        {
            Info << "iCorr ="<< iCorr << endl; 
            U.storePrevIter();
            p.storePrevIter();

            fvScalarMatrix pEqn
            (
                fvm::ddt(p)
                == 
                fvm::laplacian(Dp, p) 
                - fvc::div(fvc::ddt(Dp2,U))
            );
            pEqn.relax();
            pResidual = pEqn.solve().initialResidual();
            p.relax();

            fvVectorMatrix UEqn
            (
                fvm::laplacian(2*mu + lambda, U, "laplacian(DD,U)") 
                + divSigmaExp
                == 
                fvc::grad(p)
            );
            UEqn.relax();
            UResidual = UEqn.solve().initialResidual();
            U.relax();
            gradU = fvc::grad(U);
            sigmaD = mu*twoSymm(gradU) + (lambda*I)*tr(gradU);
            divSigmaExp = fvc::div(sigmaD - (2*mu + lambda)*gradU,"div(sigmaD)");

            residual = max(pResidual,UResidual);

        } while (residual > convergenceTolerance && ++iCorr < nCorr);

        Info << "number of iterations " << iCorr << endl;


        #include "calculateStress.H"

        #include "prepareBuildup.H"

        fvScalarMatrix pEEqn
        (
            fvm::ddt(pE) 
            == fvm::laplacian(cv, pE) 
            + f
        );
        pEEqn.solve();


        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

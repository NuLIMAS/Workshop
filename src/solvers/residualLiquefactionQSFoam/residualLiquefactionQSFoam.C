/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    residualLiquefactionQSFoam is a modified version of solidDisplacementFoam
    to predict the onset od residual liquefaction.

Description
    Transient segregated finite-volume solver of the Biot quasi-steady consolidation equations
    for linear-elastic, small-strain deformation of a solid skeleton coupled with
    pore water flow and pressure governed by Darcy's law.

    Simple linear elasticity structural analysis code.
    Solves for the displacement vector field U and the pore water pressure p.
    Also generating the stress tensor field sigma.

    An additional laplacian equation is solved in order to calculate the accumulated
    pore pressure pE.

Author
    R. Shanmugasundaram, Wikki GmbH
    H. Rusche, Wikki GmbH
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "readProperties.H"
    #include "readPoroElasticControls.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating displacement field\n" << endl;

    while (runTime.loop())
    {
        Info<< "Iteration: " << runTime.value() << nl << endl;

        #include "readPoroElasticControls.H"

        // Initialize correction iteration count and residual variables
        int iCorr = 0;
        scalar UResidual = 1.0e10;
        scalar pResidual = 1.0e10;
        scalar residual = 1.0e10;

        do
        {
            Info << "iCorr ="<< iCorr << endl;

            // Store previous iterations of U and p
            U.storePrevIter();
            p.storePrevIter();

            sigmaD.correctBoundaryConditions();

            // Solve for pressure (p) using the poroelastic equation
            fvScalarMatrix pEqn
            (
                (1/Dp2)*fvm::ddt(p)
                ==
                fvm::laplacian(Dp3, p)
                - fvc::div(fvc::ddt(U))

            );
            pEqn.relax();
            pResidual = pEqn.solve().initialResidual();
            p.relax();

            // Solve for displacement (U) using the poroelastic equation
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

        V= -Dp3 *fvc::grad(p);

        surfaceScalarField phi = - fvc::interpolate(V) & mesh.Sf();

        volScalarField contErr = fvc::div(fvc::ddt(U)) +fvc::ddt(1/Dp2,p)
                                    + fvc::div(phi);

        Info << "div1 :" <<   gSum(contErr) << endl;

        Info << "number of iterations " << iCorr << endl;

        // Include file for calculation of stress
        #include "calculateStress.H"

        // Include file for preparing buildup
        #include "prepareBuildup.H"

        // Solve the governing equation for accumulated pore pressure
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

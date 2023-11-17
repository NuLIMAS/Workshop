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
    To predict seabed liquefaction  untill  compaction

Description
    Transient segregated finite-volume solver of the Biot quasi-steady consolidation equations
    for linear-elastic, small-strain deformation of a solid skeleton coupled with
    pore water flow and pressure governed by Darcy's law.

    Simple linear elasticity structural analysis code.
    Solves for the displacement vector field U and the pore water pressure p.
    Also generating the stress tensor field sigma.

    An additional laplacian equation is solved in order to calculate the accumulated
    pore pressure pE and predicts residual liquefaction.

    Liquefied region is solved using drift-flux model.

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
    #include "readProperties.H"
    #include "readPoroElasticControls.H"
    pisoControl piso(mesh);
    #include "createFields.H"
    #include "initContinuityErrs.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    scalar sumLocalphiV = 0;

    while (runTime.loop())
    {
        Info<< "Iteration: " << runTime.value() << nl << endl;

        #include "readPoroElasticControls.H"

        // Initialize correction iteration count and residual variables
        int iCorr = 0;
        scalar UResidual = 1.0e10;
        scalar pResidual = 1.0e10;
        scalar VResidual = 1.0e10;
        scalar residual = 1.0e10;

        #include "updateValues.H"

        do
        {

            Info << "iCorr = " << iCorr << endl;
            // Store previous iterations of U, p, and V
            U.storePrevIter();
            p.storePrevIter();
            V.storePrevIter();

            // Include file for velocity equation
            #include "VEqn.H"

            // Pressure-velocity coupling loop (piso loop)
            while (piso.correct())
            {
                 // Include file for alpha transport equation (Drift-flux continuity)
                #include "alphaEqn.H"

                #include "pEqn.H"
            }

            // Include file for displacement equation
            #include "UEqn.H"

            residual = max(pResidual,UResidual);

        } while (residual > convergenceTolerance && ++iCorr < nCorr);

        Info << "Number of iterations:" << iCorr <<endl;

        // Include file for pore pressure accumulation
        #include "prepareBuildup.H"
        runTime.write();

        // Calculate and output the sum of alpha
        scalar sum = gSum(alpha);
        Info << "div1 :" <<   gSum(phiAlpha) << endl;
        std::ofstream file1;
        file1.open ("alpha.txt", std::ofstream::out | std::ofstream::app);
        file1 << runTime.timeName() << " " << sum  << std::endl << "\n";
        file1.close();

        // Include file for calculating liquefaction depth
        #include "liquefactionDepth.H"

        // Output execution time information
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
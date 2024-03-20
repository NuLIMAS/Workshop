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

    Post-liquefaction, the liquefied region is solved using drift-flux model.
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

            scalar sumLocalphiV = 0;

    while (runTime.loop())
    {
        Info<< "Iteration: " << runTime.value() << nl << endl;

        #include "updateLiquefaction.H"

        bool meshChanged = mesh.update();

        #include "readPoroElasticControls.H"
        int iCorr = 0;
        scalar UResidual = 1.0e10;
        scalar pResidual = 1.0e10;
        scalar VmResidual = 1.0e10;
        scalar residual = 1.0e10;

        #include "updateValues.H"

        do
        {
            Info << "iCorr = " << iCorr << endl;
            U.storePrevIter();
            p.storePrevIter();
            Vm.storePrevIter();

            #include "VmEqn.H"

            while (piso.correct())
            {
                #include "alphaEqn.H"

                #include "pEqn.H"
            }

            #include "UEqn.H"

            residual = max(pResidual,UResidual);

        } while (residual > convergenceTolerance && ++iCorr < nCorr);

        Info << "Number of iterations:" << iCorr <<endl;

        #include "prepareBuildup.H"

        scalar sum = gSum(alpha);
        Info << "div1 :" <<   gSum(phiAlpha) << endl;
        std::ofstream file1;
        file1.open ("alpha.txt", std::ofstream::out | std::ofstream::app);
        file1 << runTime.timeName() << " " << sum  << std::endl << "\n";
        file1.close();


        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

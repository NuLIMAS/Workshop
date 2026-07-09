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
    initialConsolidationFoam

Description
    Simple linear elasticity structural analysis code.
    Solves for the displacement vector field U, also generating the
    stress tensor field sigma.

Author
    R. Shanmugasundaram, Wikki GmbH
    H. Rusche, Wikki GmbH


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"
#   include "readMaterialProperties.H"
#   include "readSolidDisplacementFoamControls.H"
#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating displacement field\n" << endl;

// #   include "readSolidDisplacementFoamControls.H"

    // Initialize correction iteration count and residual variables
    int iCorr = 0;
    scalar UResidual = 1.0e10;
    scalar residual = 1.0e10;

    Info<< "Solving until U converges..." << nl << endl;

    do
    {
        Info << "iCorr = " << iCorr << endl;

        // Store previous iterations of U
        U.storePrevIter();

        fvVectorMatrix UEqn
        (
            fvm::laplacian(2*mu + lambda, U, "laplacian(DD,D)")
            + divSigmaExp
            //+ rhodg
        );

        UResidual = UEqn.solve().initialResidual();

        gradU = fvc::grad(U);
        strain = 0.5*twoSymm(gradU);
        volStrain = tr(strain);
        sigmaD = mu*twoSymm(gradU) + (lambda*I)*tr(gradU);

        divSigmaExp = fvc::div
        (
            sigmaD - (2*mu + lambda)*gradU,
            "div(sigmaD)"
        );
        residual = UResidual;

    } while (residual > convergenceTolerance && ++iCorr < nCorr);


    runTime++;
    # include "calculateStress.H"

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

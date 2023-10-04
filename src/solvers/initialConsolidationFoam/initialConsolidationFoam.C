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
    Solves for the displacement vector field D, also generating the
    stress tensor field sigma.
    Also, calculates the Elastic modulus and initial effective
    stress of the soil

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

    while (runTime.loop())
    {
        Info<< "Iteration: " << runTime.value() << nl << endl;

#       include "readSolidDisplacementFoamControls.H"

        int iCorr = 0;
        scalar initialResidual = 0;

        do
        {

            fvVectorMatrix UEqn
            (
                fvm::laplacian(2*mu + lambda, U, "laplacian(DD,D)")
                + divSigmaExp + rhog
            );

            initialResidual = UEqn.solve().initialResidual();

            strain = 0.5*twoSymm(gradU);
            gradU = fvc::grad(U);
            volStrain = tr(strain);
            sigmaD = mu*twoSymm(gradU) + (lambda*I)*tr(gradU);

            divSigmaExp = fvc::div
            (
                sigmaD - (2*mu + lambda)*gradU,
                "div(sigmaD)"
            );

        } while (initialResidual > convergenceTolerance && ++iCorr < nCorr);

        sigma0 = - (sigmaD.component(symmTensor::XX)
                 + sigmaD.component(symmTensor::YY)
                 + sigmaD.component(symmTensor::ZZ)) /3;

        sigmaA = sigma0 - effStress;
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

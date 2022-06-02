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
    pressureDisplacementFoam is used solved to solve pore pressure buildup
    and calculate the onset of residual liquefaction 

Description
    pressureBuildupFoam is an incremental development on biotFoam. It adds an additional
    equation (a lapalcian with source term) to solve for the pore pressure buildup. The 
    source term depends on the shear stress, initial effective stress and relative 
    density of the soil.

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
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    while (runTime.loop())
    {
        Info<< "Iteration: " << runTime.value() << nl << endl;
                            
        #include "prepareBuildup.H"
                
        fvScalarMatrix pEEqn
        (
            fvm::ddt(pE) 
            == 
            fvm::laplacian(cv, pE) + f
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

/*---------------------------------------------------------------------------*\
         _      _                 |
        | |    | |                | OLDEV: CFD tools developed by
   ___  | |  __| |  ___ __   __   |        Fraunhofer IWES and ForWind
  / _ \ | | / _` | / _ \\ \ / /   |        in Oldenburg, Germany
 | (_) || || (_| ||  __/ \ V /    |
  \___/ |_| \__,_| \___|  \_/     | http://www.iwes.fraunhofer.de
                                  | http://www.forwind.de
-------------------------------------------------------------------------------
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OLDEV and it is based on OpenFOAM.

    OLDEV is owned by the Fraunhofer IWES and ForWind Oldenburg CFD groups.
    It is not allowed to use it outside the scope of Fraunhofer IWES or
    ForWind Oldenburg CFD projects. It is forbidden to use, modify, copy or
    redistribute OLDEV files for any other purpose.

    OLDEV is developed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.

    OpenFOAM is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.


Application
    ablBoundaryFoam.C

Description

    generate the boundary condition for ablSimpleFoam in a numerical method. 
    buoyant force and Coriolis force can be switched on/off in system/fvSolution 	

Author
       Samuel Chang (chi-yao.chang@iwes.fraunhofer.de)

Modified by


Test location


Tested by


Used by


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "simpleControl.H"
#include "fvIOoptionList.H"
#include "fixedFluxPressureFvPatchScalarField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readTransportProperties.H"
    #include "createFields.H"
    #include "createGradPd.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    simpleControl simple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

        // --- Pressure-velocity SIMPLE corrector loop
        while (simple.loop())
        {
        Info<< "Time = " << runTime.timeName() << nl << endl;
         
         #include "UEqn.H"
         #include "TEqn.H"
         #include "pEqn.H"

         turbulence->correct();          



      // Correct driving force for a constant mass flow rate

        // Extract the velocity in the flow direction
        dimensionedScalar magUbarStar =
            (flowDirection & U)().weightedAverage(mesh.V());

        // Calculate the pressure gradient increment needed to
        // adjust the average flow-rate to the correct value


        dimensionedScalar gragPplus =
            (magUbar - magUbarStar)/rAU.weightedAverage(mesh.V());

        U += flowDirection*rAU*gragPplus;

        gradPd += gragPplus;

        Info<< "Uncorrected Ubar = " << magUbarStar.value() << tab
            << "correction = " <<gragPplus.value() << endl;

        runTime.write();
            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
        }



    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    writeColumns

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvc.H"
#include "IFstream.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    Foam::timeSelector::addOptions();
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createMesh.H"


    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);
  	volVectorField U_
	(
        	IOobject
        	(
        	    "U",
        	    runTime.timeName(),
        	    mesh,
        	    IOobject::MUST_READ,
        	    IOobject::NO_WRITE
        	),
        	mesh
    	);
  	volScalarField T_
	(
        	IOobject
        	(
        	    "T",
        	    runTime.timeName(),
        	    mesh,
        	    IOobject::MUST_READ,
        	    IOobject::NO_WRITE
        	),
        	mesh, 293
    	);

  	volScalarField k_
    	(
        	IOobject
        	(
        	    "k",
        	    runTime.timeName(),
        	    mesh,
        	    IOobject::READ_IF_PRESENT,
       		    IOobject::NO_WRITE
        	),
        	mesh, 0.5
    	);

    	volScalarField eps_
    	(
        	IOobject
        	(
        	    "epsilon",
        	    runTime.timeName(),
        	    mesh,
        	    IOobject::READ_IF_PRESENT,
        	    IOobject::NO_WRITE
        	),
        	mesh, 0.5
    	);


        OFstream outPutFile
        (
//            runTime.path()/"boundaryData.dat"
            "boundaryData.dat"
        );

        forAll(mesh.cells(), celli)
        {
             outPutFile	<< mesh.cellCentres()[celli][2] << tab << U_.internalField()[celli].component(0) << tab 
			<< U_.internalField()[celli].component(1) << tab << U_.internalField()[celli].component(2) << tab
		       	<< T_.internalField()[celli] << tab
			<< k_.internalField()[celli] << tab 
			<< eps_.internalField()[celli] << tab
			<< '\n';
        }
    
    }

    Info << "End\n" << endl;


    return 0;
}


// ************************************************************************* //

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

#include "Time.H"
//#include "primitiveMesh.H"
#include "argList.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{


//#   include "createMesh.H"

    
# include "setRootCase.H"
# include "createTime.H"
# include "createPolyMesh.H"
using namespace Foam;


	    IOdictionary patchNameDict
	    (
	        IOobject
	        (
	            "patchNameDict",
	            mesh.time().constant(),
	            mesh,
	            IOobject::MUST_READ_IF_MODIFIED,
	            IOobject::NO_WRITE
	        )
	    );

        const wordList patchNames(patchNameDict.lookup("patches"));

	Info<< "Dump face centres of given patch\n" << endl;
	//change the patch name to your boundary name

    forAll(patchNames, i)
    {
	label patchI = mesh.boundaryMesh().findPatchID(patchNames[i]);

        OFstream outPutFile
        (
            runTime.path()/"PatchCoordinate.dat"
        );

	forAll(mesh.boundaryMesh()[patchI].faceCentres(), faceI)
	{
		scalar x = mesh.boundaryMesh()[patchI].faceCentres()[faceI].x();
		scalar y = mesh.boundaryMesh()[patchI].faceCentres()[faceI].y();
		scalar z = mesh.boundaryMesh()[patchI].faceCentres()[faceI].z();
		outPutFile <<x<<tab<<y<<tab<<z<<tab<<endl;
	}
   }

}


// ************************************************************************* //

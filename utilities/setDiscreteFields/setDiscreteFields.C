/*---------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 1991-2010 OpenCFD Ltd.
    \\/      M anipulation   |
-------------------------------------------------------------------------------
License
    This file is based on OpenFOAM.

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
   setDiscreteFields

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "interpolateXY.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"

    instantList Times = runTime.times();
#   include "checkTimeOptions.H"
    runTime.setTime(Times[startTime],startTime);
#   include "createMesh.H"

    typedef VectorSpace<Vector<scalar>,scalar,4> quad;
    typedef VectorSpace<Vector<scalar>,scalar,6> hex;
    
    //read setDiscreteFieldsDictionary

    IOdictionary setDiscreteFieldsDict
      (
       IOobject
       (
	"setDiscreteFieldsDict",
	runTime.system(),
	mesh,
	IOobject::MUST_READ,
	IOobject::NO_WRITE
	)
       );
    // read pointer lists of each "Fields"
    PtrList<entry> parts=setDiscreteFieldsDict.lookup("Fields");

    forAll(parts,partI) {
      const dictionary &part=parts[partI].dict();
      Info <<part <<endl;

      // "field" ,"type", and "profile" are required in setDiscreteFieldsDict
      word field=part["field"];

      word type=part["type"];

      // "direction", "internal", and "patchNames" are option
      string direction="x";
      if (part.found("direction")) {
	direction=part["direction"];
      }

      bool internal=true;
      if (part.found("internal")) {
	internal=readBool(part["internal"]);
      }
      
      wordList patchNames;      
      if (part.found("patchNames")){
	patchNames=wordList(part.lookup("patchNames"));
      }

      const string &time = runTime.timeName();
      Info <<"time " <<time <<"\n";

      // for scalar field
      if(type=="scalar"){
        Info << "set " <<field <<endl;

        // read profile
	List<quad> profile(part.lookup("profile"));
	scalarField posX(profile.size());
	scalarField posY(profile.size());
	scalarField posZ(profile.size());
	scalarField fX(profile.size());
	forAll(profile,i)
	  {
	    posX[i]=profile[i][0];
	    posY[i]=profile[i][1];
	    posZ[i]=profile[i][2];
	    fX[i]=profile[i][3];
	  }
        // read target field (scalar)
	volScalarField tmpField
	  (
	   IOobject
	   (
	    field,
	    time,
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::AUTO_WRITE
	    ),
	   mesh
	   );
	
	//read cellCenter of targetMesh
	volVectorField center=tmpField.mesh().C();

	//internalField
	if(internal){
          scalarField &Ifield = tmpField.internalField();
          scalarField Icenter(Ifield.size());
          scalarField newIField(Ifield.size());
          if(direction=="x"){
            Icenter = center.internalField().component(vector::X);
            newIField = interpolateXY(Icenter,posX,fX);
          }	
          else if(direction=="y"){
            Icenter = center.internalField().component(vector::Y);
            newIField = interpolateXY(Icenter,posY,fX);
          }	
          else if(direction=="z"){
            Icenter = center.internalField().component(vector::Z);
            newIField  = interpolateXY(Icenter,posZ,fX);
          }	
          Ifield = newIField;
        }
	//patch
	forAll(patchNames,patchNameI)
	  {
	    label patchID=mesh.boundaryMesh().findPatchID(patchNames[patchNameI]);
	    scalarField &Pfield = tmpField.boundaryField()[patchID];
	    scalarField newPField(Pfield.size());
	    scalarField Pcenter(Pfield.size());
	    if(direction=="x"){
	      Pcenter = center.boundaryField()[patchID].component(vector::X);
	      newPField = interpolateXY(Pcenter,posX,fX);
	    }	
	    else if(direction=="y"){
	      Pcenter = center.boundaryField()[patchID].component(vector::Y);
	      newPField = interpolateXY(Pcenter,posY,fX);
	    }	
	    else if(direction=="z"){
	      Pcenter = center.boundaryField()[patchID].component(vector::Z);
	      newPField = interpolateXY(Pcenter,posZ,fX);
	    }	
	    Pfield = newPField;
	  }
	tmpField.write();
      }

      // for vector field
      else if(type=="vector"){
        Info << "set " <<field << endl;
        // read profile(vector)
	List<hex> profile(part.lookup("profile"));
	scalarField posX(profile.size());
	scalarField posY(profile.size());
	scalarField posZ(profile.size());
	scalarField fX(profile.size());
	scalarField fY(profile.size());
	scalarField fZ(profile.size());
	forAll(profile,i)
	  {
	    posX[i]=profile[i][0];
	    posY[i]=profile[i][1];
	    posZ[i]=profile[i][2];
	    fX[i]=profile[i][3];
	    fY[i]=profile[i][4];
	    fZ[i]=profile[i][5];
	  }

        // read target field (vector)
	volVectorField tmpField
	  (
	   IOobject
	   (
	    field,
	    time,
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::AUTO_WRITE
	    ),
	   mesh
	   );
	
	//cellCenter of targetMesh
	volVectorField center=tmpField.mesh().C();
	
	//internalField
	if(internal){
          vectorField &IVectorfield = tmpField.internalField();
          vectorField IVectorCenter(IVectorfield.size());
          vectorField newIVectorField(IVectorfield.size());
          IVectorCenter = center.internalField();
          if(direction=="x"){
            newIVectorField = interpolateXY(IVectorCenter.component(vector::X),posX,fX)*vector(1,0,0)
              +interpolateXY(IVectorCenter.component(vector::X),posX,fY)*vector(0,1,0)
              +interpolateXY(IVectorCenter.component(vector::X),posX,fZ)*vector(0,0,1);
          }	
          else if(direction=="y"){
            newIVectorField = interpolateXY(IVectorCenter.component(vector::Y),posY,fX)*vector(1,0,0)
              +interpolateXY(IVectorCenter.component(vector::Y),posY,fY)*vector(0,1,0)
              +interpolateXY(IVectorCenter.component(vector::Y),posY,fZ)*vector(0,0,1);
          }	
          else if(direction=="z"){
            newIVectorField = interpolateXY(IVectorCenter.component(vector::Z),posZ,fX)*vector(1,0,0)
              +interpolateXY(IVectorCenter.component(vector::Z),posZ,fY)*vector(0,1,0)
              +interpolateXY(IVectorCenter.component(vector::Z),posZ,fZ)*vector(0,0,1);
          }	
          IVectorfield = newIVectorField;
        }
        
	//patch
	forAll(patchNames,patchNameI)
	  {
	    label patchID=mesh.boundaryMesh().findPatchID(patchNames[patchNameI]);
	    vectorField &PVectorfield = tmpField.boundaryField()[patchID];
	    vectorField newPVectorField(PVectorfield.size());
	    vectorField PVectorCenter(PVectorfield.size());
	    PVectorCenter = center.boundaryField()[patchID];
	    if(direction=="x"){
	      newPVectorField = interpolateXY(PVectorCenter.component(vector::X),posX,fX)*vector(1,0,0)
		+interpolateXY(PVectorCenter.component(vector::X),posX,fY)*vector(0,1,0)
		+interpolateXY(PVectorCenter.component(vector::X),posX,fZ)*vector(0,0,1);
	    }	
	    else if(direction=="y"){
	      newPVectorField = interpolateXY(PVectorCenter.component(vector::Y),posY,fX)*vector(1,0,0)
		+interpolateXY(PVectorCenter.component(vector::Y),posY,fY)*vector(0,1,0)
		+interpolateXY(PVectorCenter.component(vector::Y),posY,fZ)*vector(0,0,1);
	    }	
	    else if(direction=="z"){
	      newPVectorField = interpolateXY(PVectorCenter.component(vector::Z),posZ,fX)*vector(1,0,0)
		+interpolateXY(PVectorCenter.component(vector::Z),posZ,fY)*vector(0,1,0)
		+interpolateXY(PVectorCenter.component(vector::Z),posZ,fZ)*vector(0,0,1);
	    }	
	    PVectorfield = newPVectorField;
	  }
	tmpField.write();
      }
    }
    Info << "End\n" << endl;
    return 0;
}
    
    
// ************************************************************************* //

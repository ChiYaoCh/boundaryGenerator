#! bin/sh

#function pause(){
#   read -p "$*"
#}

rm generateBoundary/boundaryData.dat inletPatchCoordinate.dat
#-------	I. generate Boundary profile
#-------I.i creating 1D mesh 		
cd generateBoundary/
echo "current directory:" $(pwd)
read -p "Press enter to continue" nothing

sh cleanResult.sh

cd ../utilities
echo "generating 1D Mesh and Profiles"
python calculateUbar.py


cd ../generateBoundary
echo "current directory:" $(pwd)
read -p "Press enter to continue" nothing

blockMesh
echo "generating boundary profiles"
echo $(pwd)
generateBoundary   
echo "writing out the boundary data"
writeBoundaryData -latestTime

cd ..
echo $pwd
echo "writing inlet patch coordinate"
writeInletPatchCoordinate

cd utilities/
echo $pwd
echo "creating the setDiscreteFieldsDict"
python writeSetDiscreteFieldsDict.py
pause 'Press [Enter] key to continue...'


cd ..
cp -r 0.ori/ 0/
setDiscreteFields


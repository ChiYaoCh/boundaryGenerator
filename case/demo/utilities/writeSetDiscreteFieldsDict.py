#! /usr/bin/python
# -*- coding: utf-8 -*-
import time, math
from scipy.optimize import fsolve
from scipy.integrate import quad
from scipy.interpolate import griddata
#from matplotlib.mlab import griddata
import numpy as np

#=====================================================================================
#
#	Author:		Samuel Chang
#	Date:		Sep. 07. 2015
#	Function:	Match the generated boundary data to 3D domain
#	Usage:		1.python writeSetDiscreteFieldsDict.py	
#
#
#*************************************************************************************

boundaryData = []
height_boundary = []
u1 = []; u2 = []; u3 = []; T = []; k = []; epsilon = [];
file= open("../generateBoundary/boundaryData.dat","r")
for line in file:
	boundaryData.append(line.split("	"))
file.close()

for line in boundaryData:
	height_boundary.append(float(line[0]))
	u1.append(float(line[1]))
	u2.append(float(line[2]))
	u3.append(float(line[3]))
	T.append(float(line[4]))
	k.append(float(line[5]))
	epsilon.append(float(line[6]))	

inletPatchXYZ = []
file= open("../inletPatchCoordinate.dat","r")
for line in file:
	inletPatchXYZ.append(line.split("	"))
file.close()
inletX = []; inletY = []; inletZ = []
for line in inletPatchXYZ:
	inletX.append(float(line[0]))
	inletY.append(float(line[1]))
	inletZ.append(float(line[2]))


u1_patch = griddata(height_boundary, np.asarray(u1), inletZ, method='linear')
u2_patch = griddata(height_boundary, np.asarray(u2), inletZ, method='linear')
u3_patch = griddata(height_boundary, np.asarray(u3), inletZ, method='linear')
T_patch = griddata(height_boundary, np.asarray(T), inletZ, method='linear')
k_patch = griddata(height_boundary, np.asarray(k), inletZ, method='linear')
epsilon_patch = griddata(height_boundary, np.asarray(epsilon), inletZ, method='linear')


#-------------	writing the interpolated data to setDiscreteFieldsDict	---------------------#


header =[]
file= open("headers/setDiscreteFieldsDictHeader","r")
for line in file:
	header.append(line)	
file.close()

headerU =[]
file= open("headers/setDiscreteFieldsDictPartU","r")
for line in file:
	headerU.append(line)	
file.close()

headerk =[]
file= open("headers/setDiscreteFieldsDictPartk","r")
for line in file:
	headerk.append(line)	
file.close()

headerT =[]
file= open("headers/setDiscreteFieldsDictPartT","r")
for line in file:
	headerT.append(line)	
file.close()

headerEps =[]
file= open("headers/setDiscreteFieldsDictPartEps","r")
for line in file:
	headerEps.append(line)	
file.close()

#-----------------writing header------------------------#
file = open("../system/setDiscreteFieldsDict","w")
for line in header:
        file.write(line)
  
#--------writing U
for line in headerU:
        file.write(line) 

for i in range(0,len(inletZ)):
	file.write("			("+str(inletX[i])+" "+str(inletY[i])+" "+str(inletZ[i])+" "+str(u1_patch[i])+" "+str(u2_patch[i])+" "+str(u3_patch[i])+")"+"\n")
file.write("                ); "+"\n")
file.write("	} "+"\n")

#--------writing k
for line in headerk:
        file.write(line) 

for i in range(0,len(inletZ)):
	file.write("			("+str(inletX[i])+" "+str(inletY[i])+" "+str(inletZ[i])+" "+str(k_patch[i])+")"+"\n")
file.write("                ); "+"\n")
file.write("	} "+"\n")

#--------writing T
for line in headerT:
        file.write(line) 

for i in range(0,len(inletZ)):
	file.write("			("+str(inletX[i])+" "+str(inletY[i])+" "+str(inletZ[i])+" "+str(T_patch[i])+")"+"\n")
file.write("                ); "+"\n")
file.write("	} "+"\n")

#--------writing epsilon
for line in headerEps:
        file.write(line) 

for i in range(0,len(inletZ)):
	file.write("			("+str(inletX[i])+" "+str(inletY[i])+" "+str(inletZ[i])+" "+str(epsilon_patch[i])+")"+"\n")
file.write("                ); "+"\n")
file.write("	} "+"\n")

file.write("); "+"\n")	

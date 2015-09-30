#! /usr/bin/python
# -*- coding: utf-8 -*-
import time, math
from scipy.optimize import fsolve
from scipy.integrate import quad
import numpy as np
#=====================================================================================
#
#	Author:		Samuel Chang
#	Date:		Jan. 29. 2015
#	Function:	generate the boundary conditions for terrain simulation using 		
#			ablBoundaryFoam and modAblSimpleFoam
#	Usage:		1.python calculateUbar.py	
#			the file (boundaryCondition) will be automatically generated
#			2. in ablBoundaryFoam folder, blockMesh it
#---------------------------------------------------------------------------------------
#	Date:		June. 18. 2015
#	updates:	1. delta for the boundary layer thickness, not for the absolute height 		
#			2. wind direction follows the regulation of meteorology rules
#	Usage:		1.python calculateUbar.py	
#			the file (boundaryCondition) will be automatically generated
#			2. in ablBoundaryFoam folder, blockMesh it
#
#
#*************************************************************************************
def Phi_m(z, Mo):
    
    if (Mo > 0):
	phi_m = -5.0*z/Mo
    elif (Mo < 0):
#	European Wind Atlas: 
#	phi_m = (1.0 - 16* z/Mo )**(0.25) -1.0

#	Panofsky
	X = (1.0 - 16*z/Mo)**0.25
	phi_m = math.log((1.0+X**2.0)/2 * ((1+X)/2.0)**2.0) -2.0* math.atan(X) + math.pi /2.0
    elif (Mo == 0):
	"""
		specific case for users convenience, Mo = 0 is here defined as neutral condition
		In general, Mo = 0 is a singular point which doesn't exist physically. 
	"""
	phi_m = 0.0		
    else :
	print "Mo is not well defined"
    
    return phi_m
#-------------------------------------------------------------------------------------
def Phi_h(z, Mo):

    if (Mo > 0):
	phi_h = -5.0*z/Mo
    elif (Mo < 0):
	phi_h = 2.0*math.log(0.5*(1.0+math.sqrt(1.0 - 16* z/Mo )))
    elif (Mo == 0):
	"""
		specific case for users convenience, Mo = 0 is here defined as neutral condition
		In general, Mo = 0 is a singular point which doesn't exist physically. 
	"""
	phi_h = 0.0
    else:
	print "Mo is not well defined"    
    
    return phi_h     
#-------------------------------------------------------------------------------------

def findUtau(Uref, Href, z0, MO):
    kappa= 0.41
    phi_m = Phi_m(Href, MO)
    utau= Uref*kappa/(math.log((Href+z0)/z0)-phi_m)	
    return utau

#-------------------------------------------------------------------------------------

def findTstar(Tref, utau, MO):

    kappa= 0.41; g =9.81
    small = 1e-10; big = 1e20 
    if (MO == 0):
	Tstar = utau*utau*Tref/big		
    else:	
	Tstar = utau*utau*Tref/(kappa*g*MO+small)	
#    kappa= 0.41
#    Tstar= (Tref-Tbot)*kappa/(math.log((Href+z0)/z0)+5.0*Href/MO)	

    return Tstar

#-------------------------------------------------------------------------------------

def calcUbar(utau, z0, zmax, zmin, MO, delta ):

    umax= calcU(zmax, utau, zmin, z0, MO, delta)
    uBar, error = quad(calcU, zmin, zmax, args=(utau, zmin, z0, MO, delta))
    uBar = uBar/(zmax-zmin)
 
    return uBar

#-------------------------------------------------------------------------------------

def calcU(z, utau,  zmin, z0, MO, delta ):

    kappa= 0.41
    if (z -zmin <= delta):
	    height = z-zmin;
    else:
            height = delta;	
    phi_m = Phi_m(height, MO)	
    u= utau/kappa*(math.log((height+z0)/z0)-phi_m)

    return u

#-------------------------------------------------------------------------------------

def calcT(z, Tref, Tstar,  zmin, z0, MO, delta ):

    kappa= 0.41
    if (z -zmin <= delta):
	    height = z-zmin;
    else:
            height = delta;
    phi_h = Phi_h(height, MO)	
    Theta =Tref+ Tstar/kappa*(math.log((height+z0)/z0)-phi_h)

    return Theta
#-------------------------------------------------------------------------------------
# warning for Mo is set as 0
def warning(MO):
	if (MO == 0):										#
		print "*******************caution**************************"
		print "you set the Monon-Obukhov number as 0"
		print "here is only the convenient way for neutral condition settlement"
		print "MO = 0 doesn't really exist"
		print "this setup is equivalent to MO -> infinity\n"
		raw_input('Press <ENTER> to continue')

#-------------------------------------------------------------------------------------
def deg2Rad(deg):
    rad = deg/180.0*np.pi
    return rad 
	
#=======================================================================================#
#				user define part					#
#=======================================================================================#
Href =100										#
Uref =10										#
windAngle = 270					# (counted from north, clockwise)	#
											#
Tref = 20+273										#	
z0 =0.006										#
											#	
Nz=80				#make sure they are identical to the 3D case		#
grad=40				#make sure they are identical to the 3D case		#				
zmax =500										#
zmin = 0.0										#
											#
delta =450										#			
MO = 0								   		#
warning(MO)										#
gCp=9.8/1006.47										#
											#
wARad= 	deg2Rad(windAngle)								#
windDir =  np.asarray([-np.sin(wARad), -np.cos(wARad), 0.0])				#										
utau= findUtau(Uref, Href, z0, MO);							#
Tstar = findTstar(Tref, utau, MO);							#
turbulenceModel = "kEpsilon"								#
#=======================================================================================#
ubar = calcUbar(utau, z0, zmax, zmin, MO, delta )
umax= calcU(zmax, utau, zmin, z0, MO, delta )
T_top= calcT(zmax, Tref, Tstar, zmin, z0, MO, delta )

print "unstable condition: MO=", MO
print "uTau:", utau
print "z0:", z0
print "******  for settlement *********"
print "uBar:", "("+str(ubar*windDir[0]), str(ubar*windDir[1]), str(ubar*windDir[2])+")"
print "u at top:", "("+str(umax*windDir[0]), str(umax*windDir[1]), str(umax*windDir[2])+")" 
print "T at bottom:", Tref
print "T at top:", T_top 

file= open("analytic.dat","w")
ubarStar=0
N=1000
file.write( 'height'+' '+'velocity(m/s)'+' '+'Temperature(K)'+'\n')
for i in range(0, N+1):
        z = zmin+(zmax-zmin)/N*i
	file.write( str(z)+' '+str(calcU(z, utau, zmin, z0, MO, delta ))+' '+str(calcT(z, Tref, Tstar, zmin, z0, MO, delta ))+'\n')
file.close()
#---------------------------------------------------IWES Owned parameters---------------------------------------------------------

if (MO > 0 and MO < 1000000 ):
	state = "stable condition"
	TEqn = "true"
#	empirical parameters for stable case
elif (MO < 0 and MO > -1000000):
	state = "unstable condition"
	TEqn = "true"
else:
	state = "neutral condition"
	TEqn = "false"
#	Turbulence and Transport Phenomena, Modelling and Simulations.	Hanjalic, 2005. p64, kEpsilon Model
#---------------------------------------------------writing the boundary file-----------------------------------------------------
header = []
file= open("headers/foamHeader","r")
for line in file:
	header.append(line)	
file.close()

file = open("boundaryConditions","w")
for line in header:
        file.write(line)  

file.write(" "+" \n") 
file.write(" "+" \n")
file.write("/*	boundary conditions for the terrain simulation  	*/"+" \n") 
file.write("/*	generated by calculateUbar.py on "+ time.strftime("%d/%m/%Y"+" at "+"%H:%M") + "	*/")
file.write(" "+" \n")
file.write(" "+" \n")
file.write("/*================="+state+"===============*/"+" \n")
file.write("/*	"+"Monin Obukhov Length:"+str(MO)+"m"+"	*/"+" \n")
file.write(" "+" \n")
file.write("/*--------------user defined settelment------------------*/"+" \n")
file.write(" "+" \n")
file.write("Href"+"	"+str(Href)+";	"+"/*Reference height*/"+" \n")
file.write("Uref"+"	"+str(Uref)+";	"+"/*Reference velocity*/"+" \n")
file.write("u_tau"+"	"+str(utau)+";	"+"/*friction velocity*/"+" \n")
file.write("delta"+"	"+str(delta)+";	"+"/*boundary layer thickness*/"+" \n")
file.write(" "+" \n")
file.write("/*-----------Boundary conditions for simulation------------------*/"+" \n")
file.write(" "+" \n")
file.write("modelName"+"	"+turbulenceModel+" \n")
file.write("z0"+"	"+"uniform "+str(z0)+";	"+"/*Roughness lenth*/"+"; \n")
file.write("zmax"+"	"+str(zmax)+";	"+"/*1D domain top position*/"+" \n")
file.write("zmin"+"	"+str(zmin)+";	"+"/*1D domain bottom position*/"+" \n")
file.write("uTop"+"	"+"("+str(umax*windDir[0])+" "+str(umax*windDir[1])+" "+str(umax*windDir[2])+")"+";	"+"/*velocity on top*/"+" \n")
file.write("uBar"+"	"+"("+str(ubar*windDir[0])+" "+str(ubar*windDir[1])+" "+str(ubar*windDir[2])+")"+";	"+"/*velocity on bottom*/"+" \n")
file.write("Tref"+"	"+str(Tref)+";	"+"/*reference temperature*/"+" \n")
file.write("Ttop"+"	"+str(T_top)+";	"+"/*pot. temperature on top*/"+" \n")
file.write("Tbot"+"	"+str(Tref)+";	"+"/*pot. temperature on top*/"+" \n")
file.write("TEqnOn"+"	"+TEqn+";	"+"/*temperature equation taken into account*/"+" \n")
file.write("Coriolis"+"	"+"false;"+"	"+"/*Coriolis force not taken into account*/"+" \n")
file.write("Omega"+"	"+"(0 0 0);"+"	"+"/*(Coriolis acceleration) (0 0 0) for no Coriolis force concerned*/"+" \n" )
file.write("Ug"+"	"+"(0 0 0);"+"	"+"/*(Geostrophic wind) (0 0 0) for no Coriolis force concerned*/"+" \n" )
file.write(" "+" \n")
file.write("/*--------------Generating blockmesh in ablBoundaryFoam------------------*/"+" \n")
file.write("Nz"+"	"+str(Nz)+";	"+"/*writing in ablBoundaryFoam/constant/polyMesh/blockMeshDict*/"+" \n")
file.write("grading"+"	"+str(grad)+";	"+"/*writing in ablBoundaryFoam/constant/polyMesh/blockMeshDict*/"+" \n")
file.write(" "+" \n")
file.write(" "+" \n")
file.write("/*	for any changes, please do it in the python code.	*/"+" \n")
file.close()




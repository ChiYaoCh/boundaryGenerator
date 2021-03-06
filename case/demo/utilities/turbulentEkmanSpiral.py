#! /usr/bin/python
# -*- coding: utf-8 -*-
import time, math
import matplotlib.pyplot as plt
import numpy as np
from numpy import arange, exp
import scipy
import scipy.stats as s
from scipy.optimize import fsolve
from scipy.integrate import quad

#=====================================================================================
#
#	Author:		Samuel Chang
#	Date:		Feb. 23. 2015
#	Function:	generate the boundary conditions for terrain simulation with Coriolis force taken into account
#			 using   		
#			
#	Usage:		1.python calculateUbar.py	
#			the file (boundaryCondition) will be automatically generated
#			2. in ablBoundaryFoam folder, blockMesh it
#
#
#*************************************************************************************
# using fsolve to define G(geostrophic wind), alpha(the angle between geostrophic to ground), gamma(the ground wind to xy system),
#			 utau(friction velocity) and zp (thickness of Prandtl layer)
def turbulentEkman(p, U, V, z, z0, f):
    G, gamma, alpha, utau, zp = p 	#alpha: transformation to arbitary angle gamma: geostrophic to ground direction changes 

    A = 1.8
    B = 4.5
    kappa=0.4
    small = 1e-10

#    print "D=", utau*zp, math.sqrt(2*kappa*utau*zp/f), "zp=", zp
    D = math.sqrt(2*kappa*utau*zp/f)
    if z <= zp :
	Ux = utau/kappa*math.log((z+z0)/z0)*math.cos(gamma)
	Vx = utau/kappa*math.log((z+z0)/z0)*math.sin(gamma)
    else:
	 
	Ux = G*(1.0-math.sqrt(2)*np.exp(-(z-zp)/D)*math.sin(gamma)*math.cos((z-zp)/D + math.pi/4.0 -gamma))   
	Vx = G*math.sqrt(2)*np.exp(-(z-zp)/D)*math.sin(gamma)*math.sin((z-zp)/D + math.pi/4.0 -gamma)    

    f1 = U -(Ux*math.cos(alpha) - Vx*math.sin(alpha))  				#  U velocity
    f2 = V -(Ux*math.sin(alpha) + Vx*math.cos(alpha)) 				#  V velocity

    f3 = G*(math.cos(gamma)-math.sin(gamma)) -utau/kappa*math.log(zp/z0)	#  Velocity on Prandtl layer height
    f4 = utau/kappa *math.sqrt((np.log(utau/(f*z0))-A)**2.0+B**2.0) - G
    f5 = zp - 0.1*D
 
    return(f1, f2, f3, f4, f5)
#-------------------------------------------------------------------------------------

def calctEkamnU(z, z0, zp, alpha, gamma, G, utau, f):

    A = -1.8
    B = 4.5
    kappa=0.4
    D = math.sqrt(2*kappa*utau*zp/f)
    if z <= zp :
	Ux = utau/kappa*math.log((z+z0)/z0)*math.cos(gamma)
	Vx = utau/kappa*math.log((z+z0)/z0)*math.sin(gamma)
    else: 
	Ux = G*(1.0-math.sqrt(2)*np.exp(-(z-zp)/D)*math.sin(gamma)*math.cos((z-zp)/D + math.pi/4.0 -gamma))   
	Vx = G*math.sqrt(2)*np.exp(-(z-zp)/D)*math.sin(gamma)*math.sin((z-zp)/D + math.pi/4.0 -gamma) 

    U = (Ux*math.cos(alpha) - Vx*math.sin(alpha))
    V = (Ux*math.sin(alpha) + Vx*math.cos(alpha))

    return U
#-------------------------------------------------------------------------------------

def calctEkamnV(z, z0, zp, alpha, gamma, G, utau, f):

    A = -1.8
    B = 4.5
    kappa=0.4
    D = math.sqrt(2*kappa*utau*zp/f)
    zp = abs(zp)
    if z <= zp :
	Ux = utau/kappa*math.log((z+z0)/z0)*math.cos(gamma)
	Vx = utau/kappa*math.log((z+z0)/z0)*math.sin(gamma)
    else: 
	Ux = G*(1.0-math.sqrt(2)*np.exp(-(z-zp)/D)*math.sin(gamma)*math.cos((z-zp)/D + math.pi/4.0 -gamma))   
	Vx = G*math.sqrt(2)*np.exp(-(z-zp)/D)*math.sin(gamma)*math.sin((z-zp)/D + math.pi/4.0 -gamma) 

    U = (Ux*math.cos(alpha) - Vx*math.sin(alpha))
    V = (Ux*math.sin(alpha) + Vx*math.cos(alpha))

    return V

#-------------------------------------------------------------------------------------

def calcUbar(zmax, zmin, z0, zp, alpha, gamma, G, utau, f):

    uBar, error = quad(calctEkamnU, 0, zmax-zmin, args=(z0, zp, alpha, gamma, G, utau, f))
    uBar = uBar/(zmax-zmin)

    return uBar
#--------------------------------------------------------------------------------------
def calcVbar(zmax, zmin, z0, zp, alpha, gamma, G, utau, f):

    	
    uBar, error = quad(calctEkamnV, 0, zmax-zmin, args=(z0, zp, alpha, gamma, G, utau, f))
    uBar = uBar/(zmax-zmin)

    return uBar
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

#=======================================================================================#
#				user define part					#
#=======================================================================================#
U = 17.0										#
V = 0.0										#
Href = 4000.0	# reference height (distance to the ground)				#
Tref = 20+273										#
latitude = 90	 # in degree (°)							#
z0 = 0.3
											#					
zmin= 266										#
zmax= Href+zmin										#
Nz=300				#make sure they are identical to the 3D case		#
grad=60				#make sure they are identical to the 3D case		#
pdType = "fixedValue"									#
Prt=1.4											#
tempEqnOn = "false"									#
CoriolisOn = "true"									#
turbulenceModel = "kEpsilonLengthLimited"						#
#=======================================================================================#
T_top=Tref
f = 1.13e-4#1.0279e-4*math.sin(latitude/180.0*math.pi)
zp0 = 60#200											
iniG = math.sqrt(U*U+V*V)
G, gamma, alpha, utau, zp =  fsolve(turbulentEkman, (iniG, math.pi/4, math.pi/4, 0.5, zp0), args=(U, V, Href, z0, f))


print "================================================================================="
print 'G:  ('+str(G*math.cos(gamma+alpha))+' '+str(G*math.sin(gamma+alpha))+' 0)'
print 'gamma: '+str(gamma*180/math.pi)+'°'+' # angle between geostrophic and ground velocity'
print 'alpha: '+str(alpha*180/math.pi)+'°'+' # angle of ground velocity'

print "neutral state with Coriolis"
print "uTau:", utau
print "z0:", z0
print "delta:", zp
	


file= open("analytic.dat","w")
for i in range(int(zmin), int(zmax)):
	z = float(i)-zmin
	ux = calctEkamnU(z, z0, zp, alpha, gamma, G, utau, f)
	uy = calctEkamnV(z, z0, zp, alpha, gamma, G, utau, f)
	file.write( '0 '+'0 '+str(i)+' '+str(ux)+' '+str(uy)+' 0 '+'\n')
file.close()

uxbar = calcUbar(zmax, zmin, z0, zp, alpha, gamma, G, utau, f)
uybar = calcVbar(zmax, zmin, z0, zp, alpha, gamma, G, utau, f)


print "******  for settlement *********"
print "uBar:", "("+str(uxbar)+' '+str(uybar)+" 0)"
print "u given top:", "("+str(U)+' '+str(V)+" 0)" 
print "utop:", "("+str(calctEkamnU(Href, z0, zp, alpha, gamma, G, utau, f))+' '+str(calctEkamnV(Href, z0, zp, alpha, gamma, G, utau, f))+" 0)"
print "================================================================================="

#---------------------------------------------------------------------------------------
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
file.write("/*	generated by turbulentEkmanSpiral.py on "+ time.strftime("%d/%m/%Y"+" at "+"%H:%M") + "	*/") 
file.write(" "+" \n")
file.write("/*=================Neutral Condition with Coriolis===============*/"+" \n")
file.write("/*	"+"Coriolis parameter:"+str(f)+"1/s"+"	*/"+" \n")
file.write(" "+" \n")
file.write(" "+" \n")
file.write("/*--------------user defined settelment------------------*/"+" \n")
file.write(" "+" \n")
file.write("Href"+"	"+str(Href)+";	"+"/*Mesoscale simulation data hight*/"+" \n")
file.write("Uref"+"	"+"("+str(U)+' '+str(V)+" 0);	"+"/*Mesoscale velocity*/"+" \n")
file.write("u_tau"+"	"+str(utau)+";	"+"/*friction velocity*/"+" \n")
file.write(" "+" \n")
file.write("/*--------------For simulation------------------*/"+" \n")
file.write(" "+" \n")
file.write("modelName"+"	"+turbulenceModel+"; \n")
file.write("z0"+"	"+"uniform "+str(z0)+";	"+"/*Roughness lenth*/"+" \n")
file.write("Omega"+"	("+str(0.0)+" "+str(f*math.cos(latitude/180.0*math.pi))+" "+str(f*math.sin(latitude/180.0*math.pi))+");	"+"/*Coriolis parameter*/"+" \n")
file.write("Ug"+"	"+'('+str(U)+' '+str(V)+' 0)'+"	"+"/*(Geostrophic wind) (0 0 0) for no Coriolis force concerned*/"+" \n" )
file.write("zmax"+"	"+str(Href)+";	"+"/*1D domain top position*/"+" \n")
file.write("zmin"+"	"+str(zmin)+";	"+"/*1D domain bottom position*/"+" \n")
file.write("uTop"+"	"+"("+str(U)+' '+str(V)+" 0);	"+"/*velocity on top*/"+" \n")
file.write("uBar"+"	"+"("+str(uxbar)+' '+str(uybar)+" 0);	"+"/*velocity on bottom*/"+" \n")
file.write("Tref"+"	"+str(Tref)+";	"+"/*reference temperature*/"+" \n")
file.write("Ttop"+"	"+str(T_top)+";	"+"/*pot. temperature on top*/"+" \n")
file.write("Tbot"+"	"+str(Tref)+";	"+"/*pot. temperature on top*/"+" \n")
file.write("Prt"+"	"+str(Prt)+";	"+"/*parameter for buoyant force*/"+" \n")
file.write("pdType"+"	"+pdType+";	"+"/*parameter for buoyant force*/"+" \n")
file.write("TEqnOn"+"	"+tempEqnOn+";	"+"/* turn on/off the temperature equation */"+" \n")
file.write("Coriolis"+"	"+CoriolisOn+";	"+"/* turn on/off the Coriolis force */"+" \n")

file.write("/*--------------Generating blockmesh in ablBoundaryFoam------------------*/"+" \n")
file.write("Nz"+"	"+str(Nz)+";	"+"/*writing in ablBoundaryFoam/constant/polyMesh/blockMeshDict*/"+" \n")
file.write("grading"+"	"+str(grad)+";	"+"/*writing in ablBoundaryFoam/constant/polyMesh/blockMeshDict*/"+" \n")
file.write(" "+" \n")
file.write(" "+" \n")

file.write(" "+" \n")

file.close()





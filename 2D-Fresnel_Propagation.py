#!/usr/bin/env python

'''
* ** *** 2D Fresnel diffraction *** ** *

python executatble script for propagation 
in Fresnel zone of data import by filename 
adress or by argument (comment or 
uncomment data = np.loadtxt(sys.argv[1])).
It assume that data are wavefront field 
describe by matrice of complexe number in 
-ascii encoding file. 
Unit in computation are micrometers (um).

author : Pierre Piault
contact : piault.pro@gmail.fr 
copyright : no copyright (2019) 
'''

import sys
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt

from math import *


#define photon feature
E = 17;                    	 # keV
lambda1 = (1.2398/E)*1e-9  	 # metres
Lambda = lambda1 *1e6		 # um
wavevector = 2*pi/Lambda;
print('wavevector (um-1) %s ' % wavevector)
print('wavelenght (um) %s ' % Lambda)

#Choose Propagation Distance 
Dmin = 1*1e6 	#um
Dmax = 21*1e6	#um
DStep = 5
D = np.linspace(Dmin,Dmax, DStep)
print(D)



#import data
data = np.loadtxt('DataFileName.dat')
#data = np.loadtxt(sys.argv[1])
I = np.abs(data)
PHI = np.angle(data)
wvf =  I * np.exp(-1j * PHI)


#define camera feature
NxCAM = wvf.shape[0]
NyCAM = wvf.shape[1]
SxCAM=SyCAM=0.63 	# um
CamXsize = NxCAM*SxCAM
CamYsize = NyCAM*SyCAM
print("camera size : %s x %s (um2)\n" %(CamXsize, CamYsize))

#padding size determination
POWx = np.log(NxCAM)//np.log(2)
POWy = np.log(NyCAM)//np.log(2)
NxPADDING = int(2**(POWx+1))
NyPADDING = int(2**(POWy+1))
print('Paddding size %s :: %s' %(NxPADDING, NyPADDING))

#make camera grid for propagator function
dfx = 1./(SxCAM*NxPADDING)
dfy = 1./(SyCAM*NyPADDING)
XX2=np.linspace(0,NxPADDING,NxPADDING, True)
XX2=XX2-np.mean(XX2)
XX2=XX2*dfx
YY2=np.linspace(0,NyPADDING,NyPADDING, True)
YY2=YY2-np.mean(YY2)
YY2=YY2*dfy
[Y2,X2]=np.meshgrid(YY2,XX2)

#data in reciprocal space
TFwvf = np.fft.fft2(wvf, (NxPADDING,NyPADDING))

#def result for each distance
Iwvf_out = np.zeros((NxCAM,NyCAM,D.size))
PHIwvf_out = np.zeros((NxCAM,NyCAM,D.size))

#Fresnel propagation forward camera at distance in array D
for (n,dval) in enumerate(D):
	print(" In progress for distance (m) : %s \n distance number : [%s:%s] \n" %(dval*1e-6,n+1,np.int(DStep)))
	
	propagator=np.exp(1j*wavevector*dval)*np.exp(-1j*pi*Lambda*dval*(X2**2 + Y2**2))
	wvf_out=np.fft.ifft2(np.fft.fftshift(TFwvf)*propagator)
	Iwvf_out[:,:,n]=np.abs(wvf_out[0:NxCAM,0:NyCAM]) 
	PHIwvf_out[:,:,n] = np.arctan2(np.real(wvf_out), np.imag(wvf_out))[0:NxCAM,0:NyCAM]
	
#show evolution with distance
for (m,dval) in enumerate (D):
	plt.figure(111)
	plt.imshow(Iwvf_out[:,:,m], cmap='gray')
	plt.title(dval*1e-6)
	plt.pause(0.1)

plt.show()

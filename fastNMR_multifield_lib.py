#Coded by D. Ban between 2013-2014
#(c) 2016 St. Jude Children's Research Hospital
#Dept. of Structural Biology

from scipy import *
from numpy import *

from scipy.optimize import leastsq
from scipy.optimize import brute
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

r_hn = 1.040E-10
theta = 17.0*(pi/180.0)
csa15N = -162.0E-6
perm_freeSPACE = 12.56637061E-7
planck = 6.62606896E-34
gammaH = 2.67513E8
gammaN = -2.7116E7

dipCONS = -1.0*(perm_freeSPACE*gammaH*gammaN*planck)/(8.0*pi*pi*r_hn*r_hn*r_hn)

def prepare_data(DATA,resNUM,resARRAY):
	out_DATA = zeros((resNUM,12))
	res_count = 0
	while res_count < resNUM:
		
		for j in range(0, len(DATA[:,0])):
		
			if DATA[j,0] == resARRAY[res_count]:
				
				for i in range(0, len(DATA[0,:])):
					out_DATA[res_count,i]= DATA[j,i]
				res_count += 1 
			if res_count == resNUM:
				break
		
	return out_DATA

def spec_simp(wX,tC,Sf):
	func = (Sf * 0.4 * tC)/(1 + (tC * wX)**2.0)
	return func

def spec_LS(wX,tC,Sf,tf):
	tau = (tC * tf)/(tC + tf)
	func = 0.4 * ((Sf * tC)/(1 + (tC * wX)**2.0) + (((1.0 - Sf)*tau)/(1 + (tau * wX)**2.0)))
	return func

def calc_nXY(specFUNC,THETA,DD,B0,delN,gN,tC,Sf,wNint):
	csaCONS = (gN*B0*delN)/sqrt(3.0)
	
	AA = 0.5 * (3.0*cos(THETA)**2.0 - 1.0)
	BB = csaCONS * DD
	
	func = (sqrt(3.0)/6.0) * AA * BB * ((4.0*specFUNC(0.0,tC,Sf)) + (3.0*specFUNC(wNint,tC,Sf)))

	return func

def calc_nXY_LS(specFUNC,THETA,DD,B0,delN,gN,tC,Sf,tf,wNint):
	csaCONS = (gN*B0*delN)/sqrt(3.0)
	
	AA = 0.5 * (3.0*cos(THETA)**2.0 - 1.0)
	BB = csaCONS * DD
	
	func = (sqrt(3.0)/6.0) * AA * BB * ((4.0*specFUNC(0.0,tC,Sf,tf)) + (3.0*specFUNC(wNint,tC,Sf,tf)))

	return func

def calc_R1N(specFUNC,DD,B0,delN,gN,tC,Sf,wNint,wHint):
	dd = DD * DD
	csaCONS = (gN*B0*delN)/sqrt(3.0)
	cc = csaCONS * csaCONS

	R1NzD = (dd/4.0) * (specFUNC(wHint-wNint,tC,Sf) + (3.0*specFUNC(wNint,tC,Sf)) + (6.0*specFUNC(wHint+wNint,tC,Sf)))

	R1NzC = cc * specFUNC(wNint,tC,Sf)

	return R1NzD + R1NzC

def calc_R1N_LS(specFUNC,DD,B0,delN,gN,tC,Sf,tf,wNint,wHint):
	dd = DD * DD
	csaCONS = (gN*B0*delN)/sqrt(3.0)
	cc = csaCONS * csaCONS

	R1NzD = (dd/4.0) * (specFUNC(wHint-wNint,tC,Sf,tf) + (3.0*specFUNC(wNint,tC,Sf,tf)) + (6.0*specFUNC(wHint+wNint,tC,Sf,tf)))

	R1NzC = cc * specFUNC(wNint,tC,Sf,tf)

	return R1NzD + R1NzC

def fit_nXYR1_ratio(params,intPARAMS,specFUNC,nXY_func,R1N_func,numRES,datARRAY):
	tauC = params[0]

	Sf = intPARAMS[0]
	THETA = intPARAMS[1]
	DD = intPARAMS[2]
	B0 = intPARAMS[3]
	delN = intPARAMS[4]
	gN = intPARAMS[5]
	wNint1000 = intPARAMS[6]
	wHint1000 = intPARAMS[7]
	wNint600 = intPARAMS[8]
	wHint600 = intPARAMS[9]
	
	LM_err = array([])
	LM_calc = array([])
	LM_data = array([])
	for j in range(0, numRES):
		errFUNC = (datARRAY[j,2]/datARRAY[j,4]) + (datARRAY[j,1]*datARRAY[j,5])/datARRAY[j,4]**2.0
		#errFUNC = 1.0
		LM_err = append(LM_err, errFUNC)

		calcNXY = nXY_func(specFUNC,THETA,DD,B0,delN,gN,tauC,Sf,wNint)
		calcR1N = R1N_func(specFUNC,DD,B0,delN,gN,tauC,Sf,wNint,wHint)

		LM_calc = append(LM_calc, calcNXY/calcR1N)
		LM_data = append(LM_data, datARRAY[j,1]/datARRAY[j,4])
	
	LSQcalc = (LM_data - LM_calc)/LM_err
	LSQvalue = sum((LM_data-LM_calc)**2.0/LM_err**2.0)
	print LSQvalue
	return LSQcalc

#def create_grid(tCrange,S2range,numRES):
#	numPARAMS = numRES + 1
#	totalCALC = tCrange*S2range*numRES
#	total_grid = zeros((totalCALC,numPARAMS))
#	cc = 0
#	while cc < totalCALC:
	
#	for i in range(0, len(tCrange)):
#		while res_count < numRES:
#			for j in range(0, len(S2range)):
#				total_grid[j,r	
		

def grid_nXY_R1N(params,intPARAMS,specFUNC,nXY_func,R1N_func,numRES,datARRAY):
	tauC = params[0]
	Sf = params[1]
	
	THETA = intPARAMS[0]
	DD = intPARAMS[1]
	B0 = intPARAMS[2]
	delN = intPARAMS[3]
	gN = intPARAMS[4]
	wNint = intPARAMS[5]
	wHint = intPARAMS[6]
	
	LM_err = array([])
	LM_calc = array([])
	LM_data = array([])
	#nXY
	for j in range(0, numRES):
		LM_err = append(LM_err, datARRAY[j,2])
		calcNXY = nXY_func(specFUNC,THETA,DD,B0,delN,gN,tauC,Sf,wNint)
		
		LM_calc = append(LM_calc, calcNXY)
		LM_data = append(LM_data, datARRAY[j,1])
	#R1N
	for j in range(0, numRES):
		LM_err = append(LM_err, datARRAY[j,5])
		calcR1N = R1N_func(specFUNC,DD,B0,delN,gN,tauC,Sf,wNint,wHint)
		LM_calc = append(LM_calc, calcR1N)
		LM_data = append(LM_data, datARRAY[j,4])
	
	LSQvalue = sum((LM_data-LM_calc)**2.0/LM_err**2.0)
#	print LSQvalue
	return LSQvalue

def grid_nXY_R1N_multi(params,intPARAMS,specFUNC,nXY_func,R1N_func,numRES,datARRAY):
	tauC = params[0]
#	Sf = params[1]

	THETA = intPARAMS[0]
	DD = intPARAMS[1]
	B0 = intPARAMS[2]
	delN = intPARAMS[3]
	gN = intPARAMS[4]
	wNint = intPARAMS[5]
	wHint = intPARAMS[6]
	
	LM_err = array([])
	LM_calc = array([])
	LM_data = array([])
	#nXY
	for j in range(0, numRES):
		Sf = params[j+1]
	#	LM_err = append(LM_err, datARRAY[j,2])
		LM_err = append(LM_err, 1.0)
		calcNXY = nXY_func(specFUNC,THETA,DD,B0,delN,gN,tauC,Sf,wNint)
		
		LM_calc = append(LM_calc, calcNXY)
		LM_data = append(LM_data, datARRAY[j,1])
	#R1N
	for j in range(0, numRES):
		Sf = params[j+1]
	#	LM_err = append(LM_err, datARRAY[j,5])
		LM_err = append(LM_err, 1.0)
		calcR1N = R1N_func(specFUNC,DD,B0,delN,gN,tauC,Sf,wNint,wHint)
		LM_calc = append(LM_calc, calcR1N)
		LM_data = append(LM_data, datARRAY[j,4])
	
	LSQvalue = sum((LM_data-LM_calc)**2.0/LM_err**2.0)
#	print LSQvalue
 #       savetxt('check_run.txt', array([LSQvalue]), fmt='%1.5f')
	return LSQvalue


def grid_nXY_R1N_multi_LS(params,intPARAMS,specFUNC,nXY_func,R1N_func,numRES,datARRAY):
	tauC = params[0]
#	Sf = params[1]
	#tf = 40.0E-12
	THETA = intPARAMS[0]
	DD = intPARAMS[1]
	B0_1000 = intPARAMS[2]
	B0_600 = intPARAMS[3]
	delN = intPARAMS[4]
	gN = intPARAMS[5]
	wNint1000 = intPARAMS[6]
	wHint1000 = intPARAMS[7]
	wNint600 = intPARAMS[8]
	wHint600 = intPARAMS[9]
	
	LM_err = array([])
	LM_calc = array([])
	LM_data = array([])
	#nXY 1000
	for j in range(0, numRES):
		Sf = params[j+1]
		tf = params[j+numRES+1]
	#	LM_err = append(LM_err, datARRAY[j,3])
		LM_err = append(LM_err, datARRAY[j,2])
		calcNXY = nXY_func(specFUNC,THETA,DD,B0_1000,delN,gN,tauC,Sf,tf,wNint1000)
		
		LM_calc = append(LM_calc, calcNXY)
		LM_data = append(LM_data, datARRAY[j,2])
	#nXY 600
	for j in range(0, numRES):
		Sf = params[j+1]
		tf = params[j+numRES+1]
	#	LM_err = append(LM_err, datARRAY[j,11])
		LM_err = append(LM_err, datARRAY[j,10])
		calcNXY = nXY_func(specFUNC,THETA,DD,B0_600,delN,gN,tauC,Sf,tf,wNint600)
		
		LM_calc = append(LM_calc, calcNXY)
		LM_data = append(LM_data, datARRAY[j,10])
		
	#R1N
	for j in range(0, numRES):
		Sf = params[j+1]
		tf = params[j+numRES+1]
	#	LM_err = append(LM_err, datARRAY[j,7])
		LM_err = append(LM_err, datARRAY[j,6])
		calcR1N = R1N_func(specFUNC,DD,B0_1000,delN,gN,tauC,Sf,tf,wNint1000,wHint1000)
		LM_calc = append(LM_calc, calcR1N)
		LM_data = append(LM_data, datARRAY[j,6])
	
	LSQvalue = sum((LM_data-LM_calc)**2.0/LM_err**2.0)
#	print LSQvalue
 #       savetxt('check_run.txt', array([LSQvalue]), fmt='%1.5f')
	return LSQvalue


#Coded by D. Ban between 2013-2014
#(c) 2016 St. Jude Children's Research Hospital
#Dept. of Structural Biology

from fastNMR_multifield_lib import *

rawDAT = genfromtxt('.csv')

wSPEC1 = 1000.3E6
wH1 = 2.0 * pi * wSPEC1
wN1 = (gammaN/gammaH)*wH1
b01 = wH1/gammaH

wSPEC2 = 600.17E6
wH2 = 2.0 * pi * wSPEC2
wN2 = (gammaN/gammaH)*wH2
b02 = wH2/gammaH
#test_nXY = calc_nXY(spec_simp,theta,dipCONS,23.5,csa15N,gammaN,48.99E-9,1.0,wN)
#test_R1N = calc_R1N(spec_simp,dipCONS,23.5,csa15N,gammaN,48.99E-9,1.0,wN,wH)

global_residues = array([])

numRES = len(global_residues)

print "num Residues"
print numRES
print "\n"

fitDAT = prepare_data(rawDAT,numRES,global_residues)
print fitDAT

p0=[208.99E-9]
setPARAMS = array([1.0,theta,dipCONS,b01,b02,csa15N,gammaN,wN1,wH1,wN2,wH2])


#fit = leastsq(fit_nXYR1_ratio,x0=p0,args=(setPARAMS,spec_simp,calc_nXY,calc_R1N,numRES,fitDAT),full_output=1)
#print fit[0]

tC_grid = arange(35,80,1)*1E-9
S2_grid = arange(0.5,1.1,0.1)

#grid_calc = zeros((len(tC_grid)*len(S2_grid)*numRES,3))
gridPARAMS = array([theta,dipCONS,b01,b02,csa15N,gammaN,wN1,wH1,wN2,wH2])


rranges = (slice(40E-9,80E-9,1.0E-9)) #tauC
range_count = 0
while range_count < len(global_residues):
	rranges = append(rranges,slice(0.6,1.0,0.01))
	range_count += 1
range_count = 0
while range_count < len(global_residues):
	rranges = append(rranges,slice(0E-12,1000E-12,100E-12))
	range_count += 1
#gridPARAMS = array([theta,dipCONS,b0,csa15N,gammaN,wN,wH])

#resbrute = brute(grid_nXY_R1N_multi, rranges, args=(gridPARAMS,spec_simp,calc_nXY,calc_R1N,numRES,fitDAT),full_output=True)
resbrute = brute(grid_nXY_R1N_multi_LS, rranges, args=(gridPARAMS,spec_LS,calc_nXY_LS,calc_R1N_LS,numRES,fitDAT),full_output=True,finish=None)
print resbrute[0]
print resbrute[1]

#lowestV = min(grid_calc[:,2])
#for c in range(0, len(grid_calc[:,0])):
#	if lowestV == grid_calc[c,2]:
#		print "SOLUTION IS "
#		print grid_calc[c]

#plt.figure(1)
#plt.plot(grid_calc[:,1],grid_calc[:,2],'ko')
#plt.show()




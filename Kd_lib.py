#Coded by D. Ban between 2013-2014
#(c) 2016 St. Jude Children's Research Hospital
#Dept. of Structural Biology

from scipy import *
from numpy import *
from scipy import linalg
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import leastsq

def create_data(residue_list, dat_pts):
	counter = 0
	output_array = zeros((len(residue_list)*dat_pts,3))
	
	while counter < len(residue_list):
		a1 = genfromtxt('res_'+repr(int(residue_list[counter]))+'_data_NPMrpl5_23oct14.txt')
		for b in range(0, len(a1[:,0])):
			where = counter*dat_pts
			output_array[where + b,0] = residue_list[counter]
			output_array[where + b,1] = a1[b,0]
			#N15_shifts
			output_array[where + b,2] = a1[b,1]
		counter += 1
	return output_array

def individual_fit(params,dat_array, current_residue,points,PT):
	Kd = params[0]
	dwMAX = params[1]
	LM_data = array([])
	LM_calc = array([])

	fit_array = zeros((points,3))	
	int_count = 0 
	for j in range(0, len(dat_array[:,0])):
		
		if current_residue == dat_array[j,0]:
			
			fit_array[int_count,0] = dat_array[j,0]
			fit_array[int_count,1] = dat_array[j,1]
			fit_array[int_count,2] = dat_array[j,2]	 
			int_count += 1

	for k in range(0, len(fit_array[:,0])):
		LM_data = append(LM_data, fit_array[k,2])
		xx = dwMAX/(2.0*PT)
		yy = (fit_array[k,1] + PT + Kd) - sqrt((fit_array[k,1] + PT + Kd)**2.0 - (4.0*fit_array[k,1]*PT))
		func = xx * yy
		LM_calc = append(LM_calc, func)

	value = (LM_data-LM_calc)
	func_value = sum((LM_data-LM_calc)**2)

	report = (params[0],params[1],func_value)

	return value 

def ind_back_plot(params,dat_array, current_residue,points,PT,xplot):
	print "HERE I AM"
	print params
	Kd = params[0]
	dwMAX = params[1]
	LM_calc = array([])

	fit_array = zeros((points,3))	
	int_count = 0 
	for j in range(0, len(dat_array[:,0])):
		
		if current_residue == dat_array[j,0]:
			
			fit_array[int_count,0] = dat_array[j,0]
			fit_array[int_count,1] = dat_array[j,1]
			fit_array[int_count,2] = dat_array[j,2]	 
			int_count += 1

	for k in range(0, len(xplot)):
		xx = dwMAX/(2.0*PT)
		yy = (xplot[k] + PT + Kd) - sqrt((xplot[k] + PT + Kd)**2.0 - (4.0*xplot[k]*PT))
		func = xx * yy
		LM_calc = append(LM_calc, func)
	    
	return LM_calc, fit_array


def global_KD(params,dat_array,PT,residue_list):
	Kd = params[0]

	LM_data = array([])
	LM_calc = array([])

	res_counter = 0
	for k in range(0, len(dat_array[:,0])):
		if residue_list[res_counter] == dat_array[k,0]:
			dwMAX = params[res_counter + 1]
		if residue_list[res_counter] != dat_array[k,0]:
			res_counter = res_counter + 1
			dwMAX = params[res_counter + 1]

		LM_data = append(LM_data, dat_array[k,2])

		xx = dwMAX/(2.0*PT)
		yy = (dat_array[k,1] + PT + Kd) - sqrt((dat_array[k,1] + PT + Kd)**2.0 - (4.0*dat_array[k,1]*PT))
		func = xx * yy
		#print "DID I DO SOMETHING?"
		LM_calc = append(LM_calc, func)

	value = (LM_data-LM_calc)
	func_value = sum((LM_data-LM_calc)**2)
    
	report = (params[0], func_value)
    
    
	#print report
	return value


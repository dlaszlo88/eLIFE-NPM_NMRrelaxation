#Coded by D. Ban between 2013-2014
#(c) 2016 St. Jude Children's Research Hospital
#Dept. of Structural Biology
from Kd_lib import *

res_list = genfromtxt('residue_list.txt')

num_of_pts = 9
Po = 258
complete_data = create_data(res_list,num_of_pts)

#print len(complete_data[0:12,:])
#p0 = [1.0, 0.1,1.2]
p0 = [0.5, 0.2]
Lplot = arange(0,800.0,1)
pp1 = PdfPages('.pdf')
res_result_ind = zeros((len(res_list), 3))
for h in range(0, len(res_list)):
	print "Res:  "+repr(res_list[h])
	fit = leastsq(individual_fit, x0=p0, args=(complete_data,res_list[h],num_of_pts,Po),full_output=1)
	print fit[0]
	res_result_ind[h,0] = res_list[h]
	res_result_ind[h,1] = fit[0][0]
	res_result_ind[h,2] = fit[0][1]

p_global = [0.01]
for f in range(0, len(res_list)):
	p_global.append(0.1)
#print p_global
#gl_fit = leastsq(global_KD, x0=p_global, args=(complete_data,Po,res_list),full_output=1)
#print gl_fit[0]
#print complete_data
gl_out_array = zeros((len(res_list),3))
for g in range(0, len(res_list)):
	gl_out_array[g,0] = res_list[g]

	plot_y, plot_points = ind_back_plot(res_result_ind[g,1:3],complete_data, res_list[g],num_of_pts,Po,Lplot)
	#plot_yGL, plot_pointsGL = ind_back_plot(gl_out_array[g,1:3], complete_data, res_list[g], num_of_pts,Po,Lplot)	
	plt.figure(g)
	plt.plot(plot_points[:,1],plot_points[:,2],'ko')
	plt.plot(Lplot,plot_y,'b-')
	#plt.plot(Lplot,plot_yGL,'r-')
	#plt.title('RES: '+repr(res_list[g])+'    Kd (uM): '+repr(round(res_result_ind[g,1],2))+' dwMAX (ppm): '+repr(round(res_result_ind[g,2],3))+'\n'+
#'KdGL (uM): '+repr(round(gl_out_array[g,1],2))+' dwMAX (ppm): '+repr(round(gl_out_array[g,2],3)))
	plt.title('RES: '+repr(res_list[g])+'    Kd (uM): '+repr(round(res_result_ind[g,1],2))+' dwMAX (ppm): '+repr(round(res_result_ind[g,2],3)))
	plt.savefig(pp1,format='pdf')
pp1.close()
print "individual fit results"
print res_result_ind
print "\n"
print "global results"
print gl_out_array

h_out = zeros((len(res_list),num_of_pts))
count_it = 0
for j in range(0, len(res_list)):
	for i in range(0, num_of_pts):
		h_out[j,i] = complete_data[count_it,2]
	
		count_it += 1

savetxt('ALLselect_titr_data.txt',h_out,fmt='%1.3f')

# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 15:28:00 2022

@author: rdiazpacheco
"""

import pandas as pd
import sys, os, argparse
import numpy as np
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters, Parameter, report_fit, Model
import peakutils
import csv as csv
from scipy.interpolate import interp1d
import xlsxwriter
import pdb # debugging tool
import re

def shotList():
	cwd = os.getcwd()
	shotlist=[]
	for subdir, dirs, files in os.walk(cwd):
		for file in files:
			path = os.path.join(subdir, file)
			if ".csv" in path:
				shotlist.append(re.findall('[0-9]{11}',path)[0])
	shotlist = sorted(shotlist, key = int) # sort shotnumbers
	return shotlist
#######
# Main #
########

def main(command_line=True):
	shot_number = int
	current_thresh = float
	baseline_min = float
	baseline_max = float
	reel_id = float
	#shot_number = int(raw_input("Input shot number:\n"))
	#reel_id = str(raw_input("Input reel_id:\n"))
	#shot_list = raw_input("Input shot list (separated by spaces): \n")
	shot_list=shotList()
	current_thresh = 10.
	baseline_min = 10
	baseline_max = 100
	plot_result = True# False
	use_calibration = False
	# Execute I-V curve
	plot_data(shot_list, baseline_min, baseline_max, current_thresh, plot_result)
    
    
#######################
# Import and fit data #
#######################
def fit_data(shot_number, baseline_min, baseline_max, current_thresh, plot_result):
	# Determine full path to the shot data file
	shot_name = str(shot_number)
	
    #data_dir = os.environ['SPARTA_DATA']
	data_file_path = shot_name + ".csv"
	
    # Extract the file header information into line strings
	header_lines = 11
	header_info = []
	f = open(data_file_path)
	for i in range(header_lines):
		header_info.append(f.readline())
		header_info[i] = header_info[i][2:-1]
	f.close()
    
	# Extract relevant header values from the strings. Note that
	# unless the user has specified the tap length on the command
	# line the tap length from the file will be used
	# (i.e. signaled by a value of tap_length<0)
    
	shot_date = header_info[0].split()[1]
	operator = header_info[1].split()[1:]
	sample_name=''
	for i in range(0,len(header_info[4].split()[1:])):
		sample_name=sample_name+header_info[4].split()[1:][i]+' '
	#tap_length = float(header_info[3].split()[4])
	tap_length = float(header_info[3].split()[1]) #<<if using new pyplate data
    
	# Load the data from the data file using pandas

	#data = pd.read_csv(data_file_path, header=5)
	data = pd.read_csv(data_file_path, header=11) #<<if using new pyplate data

	# Preserve only data above a user-specified threshold (default: 10 A)
	data_conditioned = (data[data["Shunt [A]"] > current_thresh])#<<if using new pyplate data
	#data_conditioned = (data[data["Current [A]"] > current_thresh])
    # Truncate data to calculate a user-specified baseline (default 10 A - 100 A)
	baseline_conditioned = (data_conditioned[data_conditioned["Shunt [A]"] > baseline_min])
	baseline_conditioned = (baseline_conditioned[baseline_conditioned["Shunt [A]"] < baseline_max])
	# If there is remaining data, process it
	if not data_conditioned.empty:
		# print(type(data_conditioned))
        
		# Define functional fit for V/I data
		voltages1 = data_conditioned["Tap [uV]"]
		currents1 = data_conditioned["Shunt [A]"]
		# temps1=data_conditioned["Cernox [K]"]
		# removes last data point in data (because the last data point is incorrect due to pyplate error)
		voltages = voltages1[0:-1]
		currents = currents1[0:-1]
		# temps = temps1[0:-1]
		# B=round(np.mean(data_conditioned["Hall [T]"]),2)
        
        # Calculate the voltage floor as the average of baseline array [uV]
		V_floor = np.mean(baseline_conditioned["Tap [uV]"])
        
		# Set the critical voltage [uV]
		Vc = tap_length

		# Define the fit functions
		def fit_fc(params, x, data):
			# Specify the free parameters
			Ic = params['Ic'].value
			n = params['n'].value

			# Specify the fit formula
			model = Vc*(x/Ic)**n + V_floor

		# Return the quantity that is being minimized in the fit
			return (model - data)

		def fit_fc_resistive(params, x, data):
			Ic = params['Ic'].value
			n = params['n'].value
			V_floor = params['V_floor']
			R = params['R']

			model = Vc*(x/Ic)**n + V_floor + x*R

			return model - data #that's what you want to minimize
        
        
		# Create initial guesses for the free parameters
        
		params = Parameters()
		params.add('Ic', value=36, min=1, max=1000)
		#params.add('n', value=15, min=0, max=100)
		params.add('n', value=15, min=0, max=100)
		params.add('V_floor', value= 1, min=-10, max=10)
		params.add('R', value=1e-2,min=1e-4,max=10)
        
		# Perform the fit and report on the results
		result = minimize(fit_fc_resistive, params, args=(currents, voltages))
		#result = minimize(fit_fc, params, args=(data_conditioned["Current [A]"], data_conditioned["QD Voltage [uV]"]))
		report_fit(result.params)
        
        # Get the fit result
		Ic = result.params['Ic']
		n = result.params['n']
		Ic_error = result.params['Ic'].stderr
		n_error = result.params['n'].stderr
		V_floor = result.params['V_floor']
		R = result.params['R']
		R_error = result.params['R'].stderr
		########################################## adding for plotting
		
        # Synthesize an I-V curve from the fit results
        
fit_xmin = np.min(currents)
fit_xmax = np.max(currents) * 1.05
Ic_fit_array = np.linspace(fit_xmin, fit_xmax, num=100, endpoint='true')
#V_fit_array = Vc * (Ic_fit_array / Ic) ** n + V_floor
V_fit_array = Vc*(Ic_fit_array/Ic)**n + V_floor + R * Ic_fit_array
Vc_array = np.full(100, (Vc + V_floor))
		# T=getCriticalTemp(Ic,currents,temps)
		
if plot_result:
#print"##################plotting!!!############"
			font = {'family': 'sans-serif',
					'weight': 'normal',
					'size': 15}
			plt.rc('font', **font)
​
			plot_y_min = np.min(voltages)
			plot_y_max = Vc * 2 + V_floor
​
			fig, ax = plt.subplots()
​
			plt.plot(currents, voltages, 'o', markeredgecolor='grey', markeredgewidth=1, markersize=8, mfc='none')
			plt.plot(Ic_fit_array, V_fit_array, '-', color='blue', linewidth=3)
			plt.plot(Ic_fit_array, Vc_array, '--', color='blue')
			plt.plot(np.array([Ic, Ic]), np.array([plot_y_min, plot_y_max]), '--', color='blue')
​
			ax.set_xlabel('Current (A)')
			ax.set_ylabel('Voltage (uV)')
			ax.set_title('Shot %s Ic fit' % shot_name)
​
			plt.ylim(plot_y_min, plot_y_max)
​
			# Add Tc value +/- error box to plot
			props = dict(boxstyle='round', facecolor='white', alpha=0.5)
			ax.text(0.10, 0.90, 'Sample: %s \n $I_c$ = %.2f +/- %.2f A \n n = %.2f +/- %.2f \n R=%.2f+/- %.2f nOhms' % (
				sample_name, Ic, Ic_error, n, n_error,R*1e3,R_error*1e3), transform=ax.transAxes,
					verticalalignment='top', horizontalalignment='left', bbox=props)
​
			plt.savefig(shot_name+"_re-fit.png")
​
	else:
		print("Empty dataframe!")
	return [Ic ,n ,Ic_error ,n_error , currents, voltages, V_floor, Vc, sample_name, R,R_error,shot_number]#,B,T]
​
​
######### end of fit_data
​
​
def getCriticalTemp(IC, IArray, TArray):
	setting=0.2 #K IC [A] +- 0.2 K
	start=IC-setting
	stop=IC+setting
	Temps=[]
	for i in IArray:
		if i> start and i< stop: # add Temp to list if between start and stop position
			index=IArray[IArray==i].index[0] # index = row, or index of i, element of IArray https://stackoverflow.com/questions/18327624/find-elements-index-in-pandas-series
			Temps.append(TArray[index])
	return np.mean(Temps)
​
​
​
def plot_data(shot_list, baseline_min, baseline_max, current_thresh, plot_result):
​
	print(shot_list)
	n_shots = int(len(shot_list))
	Ic_array = [0]*n_shots
	Ic_error_array = [0]*n_shots
	Vc_array = [0]*n_shots
	V_floor_array = [0]*n_shots
	n_array = [0]*n_shots
	n_error_array = [0]*n_shots	
	# B_array=[0]*n_shots
	# T_array=[0]*n_shots
	R_array=[0]*n_shots
	R_error_array=[0]*n_shots
	shot_number_array=[0]*n_shots
	currents_array = [0]*n_shots
	voltages_array = [0]*n_shots
	sample_name_array = [0]*n_shots
	Iminmax_array = np.zeros([n_shots,2])
	Vminmax_array = np.zeros([n_shots, 2])
​
	i = 0
	for shot_number in shot_list:
		try:
			# params_array = [Ic_temp, n_temp, Ic_error_temp, n_error_temp, currents_temp, voltages_temp, V_floor_temp, Vc_temp]
			params_array = fit_data(shot_number, baseline_min, baseline_max, current_thresh, plot_result)
			# print(params_array)
			Ic_array[i] = params_array[0]/1
			n_array[i] = params_array[1]/1
			Ic_error_array[i] = params_array[2]
			n_error_array[i] = params_array[3]
			currents_array[i] = params_array[4]
			voltages_array[i] = params_array[5]
			V_floor_array[i] = params_array[6]
			Vc_array[i] = params_array[7]
			sample_name_array[i] = params_array[8]
			R_array[i]=params_array[9]
			R_error_array[i]=params_array[10]
			shot_number_array[i]=params_array[11]
			# B_array[i]=params_array[9]
			# T_array[i]=params_array[10]
			Iminmax_array[i,0] = np.min(currents_array[i])
			Iminmax_array[i,1] = np.max(currents_array[i])
			Vminmax_array[i,:] = [np.min(voltages_array[i]), np.max(voltages_array[i])]
			i += 1
		except TypeError:
			pass
​
​
	# Synthesize an I-V curve from the fit results
	Iminmax_array = np.asarray(Iminmax_array)
	fit_xmin = np.min(Iminmax_array[:,0])
	fit_xmax = np.max(Iminmax_array[:,1]) * 1.05
	Ic_fit_array = np.linspace(fit_xmin,fit_xmax,num=100,endpoint='true')
	V_fit_array = np.zeros([n_shots,int(len(Ic_fit_array))]) #Vc_array*(Ic_fit_array/Ic_array)**n_array + V_floor_array
	#Vc_array2 = np.full(100,(Vc_array+V_floor_array))
	for i in range(0,n_shots):
		V_fit_array[i] = Vc_array[i]*(Ic_fit_array/Ic_array[i])**n_array[i] + V_floor_array[i]
​
	# Write to Excel
	workbk = xlsxwriter.Workbook('tape_test_xlsx.xlsx')
	worksht = workbk.add_worksheet()
	bold = workbk.add_format({'bold': True})
	# worksht.write('A1','HTS tape results', bold)
	worksht.write('A1', 'Sample Name')
	worksht.write('B1', 'Shot number')
	worksht.write('C1', 'Ic [A]')
	worksht.write('D1', 'Ic error [+/- A]')
	worksht.write('E1', 'n')
	worksht.write('F1', 'n error [+/-]')
	worksht.write('G1', 'R [nOhms]')
	worksht.write('H1', 'R error [+/-]')
​
​
	row = 1
	col = 0
	for j in range(0,n_shots):
		worksht.write(row, col, sample_name_array[j])
		worksht.write(row, col+1, shot_number_array[j])
		worksht.write(row, col+2, Ic_array[j])
		worksht.write(row, col+3, Ic_error_array[j])
		worksht.write(row, col+4, n_array[j])
		worksht.write(row, col+5, n_error_array[j])
		worksht.write(row, col+6, R_array[j]*1e3)
		worksht.write(row, col+7, R_error_array[j]*1e3)
		print("j = /n %d" % j)
		print("shot %d" % int(shot_list[j]))
		row += 1
​
	workbk.close()
​
	plot_y_min = np.min(Vminmax_array[:, 0]) - 2.
	plot_y_max = np.max(Vminmax_array[:, 1]) + 2
	plt.figure(n_shots+1)
	plt.plot(Ic_array,Vc_array,'o')
	plt.xlabel('Current [A]')
	plt.ylabel('Voltage [uV]')
​
​
###############################
# Create Python main function #
###############################
​
if __name__ == "__main__":
	main(command_line=False)
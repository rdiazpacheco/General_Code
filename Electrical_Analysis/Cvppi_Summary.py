# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 13:43:23 2023

@author: rdiazpacheco
Cavatappi Summary of Data
"""
#Tap distances
Ic_search_parameters = {
    0:["Cppi_I",200,500,4000,3830,],
    1:["Cppi_II",150,500,3400,3050]
    
    
    }

Tap_distance = {
    0:[20230207,"Cppi_I",22,15,20,70,40,15,22,20],
    1:[20230209,"Cppi_II",22,15.5,20,70,40,15.5,22,20],
    2:[20230215,"Cppi_II",22,15.5,20,70,40,15.5,22,20],
    3:[20230216,"Cppi_II",23,10,20,25,30,33,24,20],
    4:[20230227,"Cppi_I",23,10,20,25,30,33,24,20],
    5:[20230228,"Cppi_I",21,33,41,41,20,20,21,25]
    
    }


Tap_names = {
    0:[20230207,"Cppi_I-Lead (POS)", "Cppi_I-Core-1", "Cppi_I-Core-2", "Cppi_I-Core-3", "Cppi_I-Core-4", "Cppi_I-Core-5", "Cppi_I-Lead-2 (NEG)", "Cppi_I-Core-6"],
    1:[20230209,"Cppi_II-Lead (POS)", "Cppi_II-Core-1", "Cppi_II-Core-2", "Cppi_II-Core-3", "Cppi_II-Core-4", "Cppi_II-Core-5", "Cppi_II-Lead-2 (NEG)", "Cppi_II-Core-6"],
    2:[20230215,"Cppi_II-Lead (POS)", "Cppi_II-Core-1", "Cppi_II-Core-2", "Cppi_II-Core-3", "Cppi_II-Core-4", "Cppi_II-Core-5", "Cppi_II-Lead-2 (NEG)", "Cppi_II-Core-6"],
    3:[20230216,"Cppi_I-Lead (POS)", "Cppi_I-Core-1", "Cppi_I-Core-2", "Cppi_I-Core-3", "Cppi_I-Core-4", "Cppi_I-Core-5", "Cppi_I-Lead-2 (NEG)", "Cppi_I-Core-6"],
    4:[20230227,"Cppi_I-Lead (POS)", "Cppi_I-Core-1", "Cppi_I-Core-2", "Cppi_I-Core-3", "Cppi_I-Core-4", "Cppi_I-Core-5", "Cppi_I-Lead-2 (NEG)", "Cppi_I-Core-6"],
    5:[20230228,"Cppi_I-Lead (POS)", "Cppi_I-Core-1", "Cppi_I-Lead-1.1 (POS)", "Cppi_I-Lead-1.2 (POS)", "Cppi_I-Core-4", "Cppi_I-Core-5", "Cppi_I-Lead-2 (NEG)", "Cppi_I-Core-6"]
    }
#"""

ramp_rates_ls = {
    0:[20230207,"Cppi_I",50,50,50,50,50,50],
    1:[20230209,"Cppi_II",30,30,30,30,50,50,50,50,50,75,75,75],
    2:[20230215,"Cppi_II",50,30,50,100,100,100,30,30,200,200],
    3:[20230216,"Cppi_II",30,30,50,50,100,100,200,200],
    4:[20230227,"Cppi_I",50,50,100,100,200,200],
    5:[20230228,"Cppi_I",50,50,100,100,200,200]
    }

ramp_steps = {
    0:[20230207,"Cppi_I",3500,3700,3900,4000,4250],
    1:[20230209,"Cppi_II",0],
    2:[20230215,"Cppi_II",2700,2800,2900,3000,3250],
    3:[20230216,"Cppi_II",2700,2800,2900,3000,3250],
    4:[20230227,"Cppi_I",150,3700,3800,3900,4000],
    5:[20230228,"Cppi_I",150,3700,3800,3900,4000]
    }

ramp_rates_steps = {
    0:[20230207,"Cppi_I",50,50,50],
    1:[20230209,"Cppi_II",0],
    2:[20230215,"Cppi_II",30,50,100,200],
    3:[20230216,"Cppi_II",200,200,100,50,30],
    4:[20230227,"Cppi_I",200,200,100,50,100],
    5:[20230228,"Cppi_I",200,200,100,50]
    }

# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 15:04:34 2023

@author: rdiazpacheco


SULTAN Summary
"""



Tap_name = {
    1:["POS Lead-A Ind Long","SULTAN PVJ S-1", "SULTAN PVJ 1","SULTAN PVJ 1"],
    2:["POS Lead-B Ind Long","SULTAN PVJ S-4","SULTAN PVJ 2","SULTAN PVJ 2"],
    3:["SULTAN PVJ-A","SULTAN PVJ-A","SULTAN PVJ 3","SULTAN PVJ 3"],
    4:["SULTAN PVJ-B","SULTAN PVJ-B","SULTAN PVJ 4","SULTAN PVJ S1"],
    5:["SULTAN PVJ-C","SULTAN PVJ-C", "SULTAN PVJ 5","SULTAN PVJ S2"],
    6:["SULTAN Overall","SULTAN PVJ S-2", "SULTAN Indigo Long","SULTAN PVJ S3"],
    7:["NEG Lead-A Ind Short","SULTAN PVJ S-3", "SULTAN Lead B","SULTAN Lead C"],
    8:["NEG Lead-B Ind Short","Indigo Short", "SULTAN Lead A", "SULTAN Lead D"]
    }
#"""


Tap_dist = {
    1:[512,305,215,215],
    2:[512,455,282.5,282.5],
    3:[220,220,360,360],
    4:[280,280,422.5,1197],
    5:[520,520,1750,1217],
    6:[1930,355,370,1237],
    7:[512,405,642.5,508],
    8:[512,100,642.5,508]
    }
    

Tap_color = {
    1: "tab:blue",
    2: "tab:red",
    3: "tab:red",
    4: "tab:orange",
    5: "tab:pink",
    6: "tab:blue",
    7: "tab:blue",
    8: "tab:red"
    }

ramp_steps = {
    0:[200, 1000, 1500, 2000, 2500, 3000, 3500, 4000]
    }

ramp_rates_steps = {
    1:100,
    2:100,
    3:100,
    4:100,
    5:100
    }

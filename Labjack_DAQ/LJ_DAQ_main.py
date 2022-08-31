# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 07:56:02 2022

@author: rdiazpacheco

Based on Dan Nash Labjack Program author: dnash
"""

import sys
import json
import queue
import time

from tkinter import filedialog
from labjack_sampler import Sampler
from influx_logger import Logger
#from cryocon_sampler import CC_Sampler

if __name__ == '__main__':
      
    # configure InfluxDB connection, database, and measurement information
    with open(filedialog.askopenfilename(title = "Select InfluxDB Configuration File...",filetypes = (("JSON files","*.json"),("all files","*.*")))) as f:
        influx_config = json.load(f)
        print('Imported InfluxDB configuration from %s.' % f.name)
    
    # open JSON config file for LabJack, includes device and channel information
    with open(filedialog.askopenfilename(title = "Select Labjack Configuration File...",filetypes = (("JSON files","*.json"),("all files","*.*")))) as f:
        labjack_config = json.load(f)
        print('Imported LabJack configuration from %s.' % f.name)
    
    # determine if a cryocon 24C will be sampled from
    
    # interthread communication between Sampler() and Logger()
    data_pipe = queue.Queue()
    
    # instantiate threads
    sampler = Sampler(config=labjack_config, dout=data_pipe, verbose=True)
    logger = Logger(config=influx_config, din=data_pipe, verbose=True)
    
    # start threads
    sampler.start()
    logger.start()
    
    
    # run until KeyboardInterrupt
    try:
        while True:
            pass
    except KeyboardInterrupt:
        # trigger threads to internally end their run() loops
        sampler.kill()
        logger.kill()
     
        
        time.sleep(2)
        
        sampler.join()
        logger.join()
        
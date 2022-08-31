# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 11:48:31 2022

@author: rdiazpacheco
Based on Dan Nash
"""


import time
import threading
import queue
import json
import copy

import numpy as np

from labjack import ljm

class Sampler(threading.Thread):
    '''The Sampler() class provides the main interface to streaming JabJack data.
    '''
    def __init__(self, config=None, verbose=False, dout=None):

        super(Sampler, self).__init__()

        # LabJack T7 connection status
        self.connected = False

        # thread stuff
        self.t_name = 'SAMP'    # thread name to display in CLI
        self.terminate = False  # raise flag which ends run()
        self.verbose = verbose  # for verbose printing

        # create queue for sending data to influx logger thread
        if dout == None:
            dout = queue.Queue()
        self.dout = dout

        # Sample rate stuff
        self._fs = config['fs']                     # sampling frequency [Hz]
        self._dt = 1./self._fs

        print('%s:\tCreated sampler with %i Hz sample rate' % (self.t_name, self._fs))

        # initialize and connect to LabJack
        self.labjack = self.init_labjack(config)
        self.channels = config['channels']          # dictionary of channel information
        self.scan_list = list(self.channels.keys()) # list of LabJack / MODBUS-recognized channel names

    def init_labjack(self, config, handle=None):
        '''Connect to and configure the LabJack using device and channel information from labjack_config.
        '''
        try:
            handle = ljm.openS(deviceType=config['deviceType'], \
                               connectionType=config['connectionType'], \
                               identifier=config['identifier'])
            info = ljm.getHandleInfo(handle)
            print('%s:\tOpened a LabJack with Device type: %i, Connection type: %i, Serial number: %i, IP address: %s, Port: %i, Max bytes per MB: %i' % \
                  (self.t_name, info[0], info[1], info[2], ljm.numberToIP(info[3]), info[4], info[5]))
            self.connected = True

            # configure channels
            self.scan_addrs = []    # blank list for storing byte offsets for each channel
            settings = ['NEGATIVE_CH', 'RANGE', 'RESOLUTION_INDEX', 'SETTLING_US']     # settings by name
            for channel in config['channels']:      # each key is a LJM 'name' e.g. 'AIN0'
                values = []     # values of each setting name
                names = []      # name of each setting
                addr = config['channels'][channel]['index']
                self.scan_addrs.append(2*addr)      # byte offsets for each index
                for setting in settings:
                    name = '%s_%s' % (channel, setting)             # concatenate channel and setting name, e.g. 'AIN0_RANGE'
                    val = config['channels'][channel][setting]      # look up value from config dictionary e.g. {'AIN0': {'RANGE': 0.0'}}

                    # build list of name-value pairs to write with eWriteNames()
                    names.append(name)
                    values.append(val)
                    if self.verbose:
                        print(name,val)
                ljm.eWriteNames(handle, len(names), names, values)  # configure each channel

        except TypeError as e:
            print('%s:\tCould not connect to LabJack, with error: %s' % (self.t_name, e))

        return handle

    def pack(self, data, scan_rate, t0, i_next=0):
        '''Read in a full scan from eStreamRead() and parse into individual samples by channel. Format and put into dout queue
        
        data: list of length n_chans * scans_per_read, ordered by channel and sample.
        
        t0: time that stream started, from which to back out relative sample time.
        '''
        # template dictionary to clone for each sample
        subframe = {'time': None,
                    'fields': {}}
        outs = []   # list of samples to write to influx at one time, cleared after write

        n_chans = len(self.scan_list)
        i_samp = i_next # sample indexer, start at latest index and increment to i_last + (len(data)/n_chans)
        while data:     # keep popping until data list is empty
            '''Crawl through each individual sample, determine index, apply calibration.
            '''
            dt = (1/scan_rate)*i_samp*1e9   # convert relative time to [ns]
            t = int(t0 + dt)
            out = copy.deepcopy(subframe)
            out['time'] = t
            for i in range(n_chans):
                mod_name = self.scan_list[i]                # lookup MODBUS channel name, e.g. 'AIN0'
                log_name = self.channels[mod_name]['name']  # lookup signal name to log as, e.g. 'pressure_0'
                m = self.channels[mod_name]['slope']        # lookup calibration slope
                b = self.channels[mod_name]['offset']       # lookup calibration offset
                x = data.pop(0)                             # pop the 0th index out of the data FIFO
                out['fields'][log_name]= m*x+b              # assign a scaled value to log_name and put in the outbound dict
            outs.append(out)
            i_samp += 1

        # dump list of outbound dicts into dout queue to be picked up by influx logger
        self.dout.put(outs)
        # pass the index of the latest sample to be used for the next data FIFO to maintain timing awareness
        return i_samp

    def kill(self):
        '''Send trigger to run() to exit main loop.
        '''
        print('%s:\tStopping thread...' % self.t_name)
        self.terminate=True

    @property
    def fs(self):
        '''Sampling frequency.
        '''
        return self._fs

    @fs.setter
    def fs(self, fs):
        dt = 1./fs
        self._dt = dt
        self._fs = fs

    @property
    def dt(self):
        '''Sample time.
        '''
        return self._dt

    def run(self, mode='stream'):
        '''Sampler thread __main__ routine.
        
        If mode == 'stream', initialize and start streaming with ljm.eStreamStart(), and read from LJM buffer with ljm.eStreamRead().
        '''
        if mode == 'stream':
            # configure eStreamRead()
            n_chans = len(self.scan_list)
            scan_addrs = self.scan_addrs
            scan_rate = self._fs                # scan frequency
            scans_per_read = int(scan_rate/4)   # scans per read (0.25s)
    
            # ensure triggered stream is disabled.
            ljm.eWriteName(self.labjack, "STREAM_TRIGGER_INDEX", 0)
            # enabling internally-clocked stream.
            ljm.eWriteName(self.labjack, "STREAM_CLOCK_SOURCE", 0)
            # begin streaming with eStreamStart()
            ljm.eStreamStart(self.labjack, scans_per_read, n_chans, scan_addrs, scan_rate)

            # initialize array of correct length to be populated by eStreamRead()
            ret = np.array(n_chans*scans_per_read, dtype=float)
            tik = time.time_ns()        # record stream start time in [ns]

            i_next = 0      # most recent sample read and logged
            while not self.terminate:
                '''Main streaming routine using eStreamRead().
                '''
                t = time.time_ns()-tik      # record relative time in [ns]
                ret = ljm.eStreamRead(self.labjack)
                data = ret[0]
                samps = len(data)/n_chans
                i_last = self.pack(data, scan_rate, tik, i_next)
                if self.verbose:
                    print('%s:\tScan %s start %0.3fs\tScans: %i\tScan Backlogs: Device: %i, LJM: %i' % (self.t_name, i_next, t/1e9, samps, ret[1], ret[2]))
                i_next = i_last

        ljm.close(self.labjack)
        print('%s:\tThread successfully ended.' % self.t_name)

if __name__ == '__main__':
    # open JSON config file for LabJack
    with open('config_files/labjack_config.json','r') as f:
        labjack_config = json.load(f)
    sampler = Sampler(config=labjack_config, verbose=True)

    sampler.start()

    try:
        while True:
            pass
    except KeyboardInterrupt:
        sampler.kill()
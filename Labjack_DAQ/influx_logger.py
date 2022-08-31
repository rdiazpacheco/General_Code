# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 11:50:05 2022

@author: rdiazpacheco
From Dan Nash
"""




import time
import threading
import queue

from influxdb import InfluxDBClient

class Logger(threading.Thread):
    
    def __init__(self, config=None, din=None, verbose=False):
        
        super(Logger, self).__init__()
        self.t_name = 'LOG'

        self.config = config
        self.client = InfluxDBClient(host=config['HOST'], \
                                     port=config['PORT'], \
                                     username='root', \
                                     password='root', \
                                     database=config['DB'])
        self.meas = config['MEAS']

        # queue for receiving incoming 
        if din == None:
            din = queue.Queue()
        self.din = din
        self.terminate = False

        self.verbose = verbose

    def kill(self):
        '''Send trigger to run() to exit main loop.
        '''
        print('%s:\tStopping thread...' % self.t_name)
        self.terminate=True

    def run(self):
        '''Logger thread __main__ routine.
        '''
        while not self.terminate:
            if not self.din.empty():        # check din for available incoming messages to write to influx
                outs = self.din.get()
                for samp in outs:
                    samp['measurement'] = self.meas     # dictionaries coming from din have 'time' and 'fields', but no 'measurement' or 'tags'
                if self.verbose:
                    print('%s:\tWriting new sample to %s.' % (self.t_name, self.config['DB']))
                self.client.write_points(outs, time_precision='n')
            time.sleep(0.02)
        print('%s:\tThread successfully ended.' % self.t_name)

if __name__ == '__main__':
    influx_config = {'HOST': 'blazar',
                     'PORT': 8086,
                     'DB': 'sparc_data',
                     'MEAS': 'test'}
    
    logger = Logger(influx_config)
    
    logger.start()
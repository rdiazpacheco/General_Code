# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 11:49:06 2022

@author: rdiazpacheco
From Dan Nash
"""
# -*- coding: utf-8 -*-
import json

filename = 'labjack_config.json'

labjack_config = {'deviceType': 'T7',
                  'connectionType': 'ETHERNET',
                  'identifier': '172.20.200.25',
                  'channels': {'AIN0': {'name': 'test_se_0',
                                        'index':0,
                                        'units': 'volts',
                                        'slope': 1.0,
                                        'offset': 0.0,
                                        'NEGATIVE_CH': 199,
                                        'RANGE': 0,
                                        'RESOLUTION_INDEX': 0,
                                        'SETTLING_US': 0},
                               'AIN1': {'name': 'test_se_1',
                                        'index':1,
                                        'units': 'volts',
                                        'slope': 1.0,
                                        'offset': 0.0,
                                        'NEGATIVE_CH': 199,
                                        'RANGE': 0,
                                        'RESOLUTION_INDEX': 0,
                                        'SETTLING_US': 0},
                               'AIN2': {'name': 'test_diff_2',
                                        'index':2,
                                        'units': 'volts',
                                        'slope': 1.0,
                                        'offset': 0.0,
                                        'NEGATIVE_CH': 3,
                                        'RANGE': 0,
                                        'RESOLUTION_INDEX': 0,
                                        'SETTLING_US': 0},
                               'AIN3': {'name': 'test_se_3',
                                        'index':3,
                                        'units': 'volts',
                                        'slope': 1.0,
                                        'offset': 0.0,
                                        'NEGATIVE_CH': 199,
                                        'RANGE': 0,
                                        'RESOLUTION_INDEX': 0,
                                        'SETTLING_US': 0}}}

config_json = json.dumps(labjack_config, indent=4)

with open(filename, 'w+') as f:
    f.write(config_json)

print(config_json)
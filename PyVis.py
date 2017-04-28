# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 16:20:46 2017

@author: admin
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math, os

cwd = os.getcwd();
dataFolder = os.path.join(cwd,'data')
dataFile = os.path.join(dataFolder,'obj_pose-laser-radar-synthetic-input.txt')
outFolder = os.path.join(cwd, 'UKF')
outFile = os.path.join(outFolder,'output.txt')

#out_cols = ['time_stamp','px_state','py_state','v_state','yaw_angle_state',
#           'yaw_rate_state','sensor_type','NIS','px_measured','py_measured',
#           'px_ground_truth','py_ground_truth','vx_ground_truth','vy_ground_truth']
with open(outFile) as f:
    #data = pd.read_table(f, sep='\t', header=None, names=out_cols, lineterminator='\n')
    data = pd.read_table(f, sep='\t', lineterminator='\n')
    
# Start with plotting

# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 16:20:46 2017

@author: admin
"""

import pandas as pd
import numpy as np
from scipy.stats import chi2
import matplotlib.pyplot as plt
import os

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
plt.figure(figsize=(20,20))
ax = plt.subplot(2,2,1)
ax.set_title('XY Position')
p_gt = ax.scatter(data['px_ground_truth'], data['py_ground_truth'], color='orange', label = 'p_ground_truth')
p_state = ax.scatter(data['px_state'], data['py_state'], color='b', label = 'p_state')
#ax.scatter(data['px_measured'], data['py_measured'], color = 'r')
ax.legend(handles = [p_gt,p_state])

vx = data['v_state'] * np.cos(data['yaw_angle_state'])
vy = data['v_state'] * np.sin(data['yaw_angle_state'])

ax = plt.subplot(2,2,2)
ax.set_title('XY Position')
v_gt = ax.scatter(data['vx_ground_truth'], data['vy_ground_truth'], color='orange', label = 'v_ground_truth')
vel = ax.scatter(vx, vy, color='b', label='v_state')
ax.legend(handles=[v_gt,vel])

# Getting NIS data and the chi^2 distribution value
nis_radar = data.NIS[data.sensor_type=='radar'].as_matrix()
nis_lidar = data.NIS[data.sensor_type=='lidar'].as_matrix()

chi_radar = np.full(nis_radar.shape, chi2.isf(q=0.05, df=3))
chi_lidar = np.full(nis_lidar.shape, chi2.isf(q=0.05, df=2))

ax = plt.subplot(2,2,3)
ax.set_title('RADAR NIS')
ax.plot(nis_radar,'b')
ax.plot(chi_radar,'orange')

ax = plt.subplot(2,2,4)
ax.set_title('LIDAR NIS')
ax.plot(nis_lidar,'b')
ax.plot(chi_lidar,'orange')

# For debugging
#ax = plt.subplot(2,3,5)
#ax.set_title('Phi')
#ax.plot(data.yaw_angle_state,'b')

#ax = plt.subplot(2,3,6)
#ax.set_title('Phi_dot')
#ax.plot(data.yaw_rate_state,'b')

plt.tight_layout()
#plt.show()
plt.savefig('post_processing.jpg')
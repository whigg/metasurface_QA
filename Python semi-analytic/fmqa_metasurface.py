# -*- coding: utf-8 -*-
"""
Created on Mon May 11 15:13:02 2020

@author: cr26
"""




import numpy as np
import dimod
import fmqa

Eo = 1; 
k_x = 1;
k_y = 1;
d=0.049;
c=3e8;
    #Number of elements in x-y directions in reflectarray
 
wavelength = 0.049
pi = np.pi
k = 2*pi/wavelength

N = 6 
A=(N*d)**2
#opt_code = [[0,0,1,1,0,1,0]]  
sample_size = 10

xs = np.array(np.random.randint(2, size=(sample_size,N)))
ys = np.array([metasurface_DBS(N,x) for x in xs])

model = fmqa.FMBQM.from_data(xs, ys)

sa_sampler = dimod.samplers.SimulatedAnnealingSampler()

for _ in range(30):
    res = sa_sampler.sample(model, num_reads=3)
    xs = np.r_[xs, res.record['sample']]
    ys = np.r_[ys, [metasurface_DBS(N,x) for x in res.record['sample']]]
    model.train(xs, ys)
    
import matplotlib.pyplot as plt
plt.plot(ys, 'o')
plt.xlabel('Selection ID')
plt.ylabel('value (scaled)')
plt.show()
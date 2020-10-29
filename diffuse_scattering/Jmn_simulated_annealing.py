# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 07:47:35 2020

@author: cr26
"""
import dimod
import neal

def Jmn_int(m1,n1,m2,n2):
    import numpy as np
    from scipy import integrate
    wavelength = 0.029979
    pi = np.pi
    k = 2*pi/wavelength
    kx = lambda t,p: k*np.sin(t)*np.cos(p)
    ky = lambda t,p: k*np.sin(t)*np.sin(p)
    d=0.049;
    P = lambda t,p: np.cos(t)**2*np.sin(p)**2+np.cos(p)**2*np.sinc(k*d/2*np.sin(t)*np.cos(p))**2*np.sinc(k*d/2*np.sin(t)*np.sin(p))**2 
    integrand = lambda t,p: np.cos((m2-m1)*kx(t,p)*d+(n2-n1)*ky(t,p)*d)*P(t,p)*np.sin(t)
    J = integrate.dblquad(integrand,0,pi/2,0,2*pi)
    return J[0]

def Jmn_sum(m1,n1,m2,n2):
    theta = np.linspace(0,np.pi/2,100)
    phi = np.linspace(0,2*np.pi,100)
    wavelength = 0.029979
    k = 2*np.pi/wavelength
    kx = lambda t,p: k*np.sin(t)*np.cos(p)
    ky = lambda t,p: k*np.sin(t)*np.sin(p)
    d=2*wavelength
    J = 0
    P = lambda t,p: np.cos(t)**2*np.sin(p)**2+np.cos(p)**2*np.sinc(k*d/2*np.sin(t)*np.cos(p))**2*np.sinc(k*d/2*np.sin(t)*np.sin(p))**2 
    summand = lambda t,p: np.cos((m2-m1)*kx(t,p)*d+(n2-n1)*ky(t,p)*d)*P(t,p)*np.sin(t)
    for t in theta:
        for p in phi:
            J+=summand(t,p)
    return J



N=6
h={}
for i in range(N):
    m_str = "x"+str(i)
    h[m_str] = 0

J = {}
for m1 in range(0,N):
    for n1 in range(0,N):
        for m2 in range(0,N):
            for n2 in range(0,N):
                if(m1 != m2 or n1 != n2):
                    str_mn = ('x'+str(m1*N+n1),'x'+str(m2*N+n2))
                    J[str_mn]  = (Jmn_int(m1,n1,m2,n2))

import numpy as np
sampleset = dimod.SimulatedAnnealingSampler().sample_ising(h, J)
print(sampleset.first)
mat_dict = sampleset.first[0]
opt_code_2d = np.zeros(shape=(N,N))
for m in range(N):
    for n in range(N):
        dict_string = 'x'+str(m*N+n)
        opt_code_2d[m,n] = mat_dict[dict_string]

RCS_calc(opt_code_2d,N)
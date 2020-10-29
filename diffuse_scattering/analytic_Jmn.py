# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 10:09:54 2020

@author: cr26
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

def R(N,n):
    phase = np.zeros(N)
    phase[n] = 1
    phase = phase.reshape(int(N**.5),int(N**.5))
    phase = np.repeat(phase,10,axis=0)
    phase = np.repeat(phase,10,axis=1)
    R_freq = (np.fft.fft2(phase))
    return R_freq
def Jmn(N,m,n):
    if m==n:
        return 0
    else:
        np.exp(-1j*(k*d*np.sin(t)*((m-0.5)*np.cos(p)+(n-0.5)*np.sin(p))))
        dot_product = R(N,m)
        return abs(sum(sum(dot_product)))
def Jmn_int(m1,n1,m2,n2):
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

def Jmn_int_1d(m,n):
    wavelength = 0.029979
    pi = np.pi
    k = 2*pi/wavelength
    kx = lambda t,p: k*np.sin(t)*np.cos(p)
    ky = lambda t,p: k*np.sin(t)*np.sin(p)
    d=0.049;
    P = lambda t,p: np.cos(t)**2*np.sin(p)**2+np.cos(p)**2*np.sinc(k*d/2*np.sin(t)*np.cos(p))**2*np.sinc(k*d/2*np.sin(t)*np.sin(p))**2 
    integrand = lambda t,p: np.cos((m-n)*kx*d)*P(t,p)*np.sin(t)
    J = integrate.dblquad(integrand,0,pi/2,0,2*pi)
    return J[0]
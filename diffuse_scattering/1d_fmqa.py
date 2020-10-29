# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 11:38:57 2020

@author: cr26
"""# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 21:05:56 2020

@author: cr26
"""
def array_factor(N,kx,ky,gamma):
    arr_factor = 0
    d=0.049;
    for x in range(N):
        for y in range(N):
            arr_factor+= gamma[x,y] * np.exp(-1j*d*(kx*x+ky*y))
    return arr_factor

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

def metasurface_DBS1(N,opt_code):

    Eo = 1; 
    k_x = 1;
    k_y = 1;
    d=0.049;
    c=3e8;
        #Number of elements in x-y directions in reflectarray
     
    wavelength = 0.029979
    pi = np.pi
    k = 2*pi/wavelength

    A=(N*d)**2
    
    Es_angular_y = np.zeros(shape=(100,100),dtype=np.complex_)
    Es_angular_z = np.zeros(shape=(100,100),dtype=np.complex_)    
    D = np.zeros(shape=(100,100),dtype=np.complex_)   
    f = np.zeros(shape=(100,100),dtype=np.complex_) 
    F = np.zeros(shape=(100,100),dtype=np.complex_) 
    E_freq = np.zeros(shape=(N,N),dtype=np.complex_)
    phase = np.zeros((N,N))
    
    matrix = np.repeat([opt_code],N,axis=0)
    for i in range(N):
        if(matrix[0,i] == 1):
            matrix[i] = ~matrix[i] +2
    
    for i in range(N):
        for j in range(N):
            phase[i,j] = matrix[i,j]
    gamma = 2*phase-1
    
    #plt.imshow(phase)
    #plt.show()
    Es_freq = np.zeros((N,N))
    #for x in range(N):
    #    for y in range(N):
    #        E_temp = np.zeros(shape=(N,N),dtype=np.complex_)
    #        E_temp[x,y] = phase[10*x:10*x+10,10*y:10*y+10]*gamma[10*x:10*x+10,10*y:10*y+10]
    #        Es_freq = Es_freq + np.fft.fftshift((np.fft.fft2(E_temp)))
    
    #plt.imshow(abs((Es_freq)))
    #plt.show()
    
    Esy = lambda t,p: Eo*d**2*np.sinc(k*np.sin(t)*np.cos(p)*d/2)*np.sinc(k*np.sin(t)*np.sin(p)*d/2)
    Es_theta_y = lambda t, p: Esy(t,p)*(np.cos(t))
    Es_theta_z = lambda t,p: Esy(t,p)*(np.sin(t)*np.sin(p))
    Es_theta_squared = lambda t,p: abs(np.sqrt(Es_theta_y(t,p)**2+Es_theta_z(t,p)**2))**2
    D_denominator = lambda t,p : np.sin(t)*Es_theta_squared(t,p)*abs(array_factor(N,k*np.sin(t)*np.cos(p),k*np.sin(t)*np.sin(p),gamma))**2
    D_numerator = lambda t,p: 4*pi*Es_theta_squared(t,p)*abs(array_factor(N,k*np.sin(t)*np.cos(p),k*np.sin(t)*np.sin(p),gamma))**2
    
    #D_den = integrate.dblquad(D_denominator,0,2*pi,0,pi/2,epsabs=1e-5)
    theta = np.linspace(0,pi/2,100)
    phi = np.linspace(0,2*pi,100)
    i=0
    j=0
    for t in theta:
        for p in phi:
            f[i,j] = Es_theta_squared(t,p)
            Es_angular_y[i,j] = Es_theta_y(t,p)
            Es_angular_z[i,j] = Es_theta_z(t,p)
            F[i,j] = array_factor(N,t,p,gamma)
            D[i,j] = D_numerator(t,p)#/D_den[0]
            i+=1
        j+=1
        i=0
        
    #plt.imshow(D.real)
    RCS = wavelength**2/(4*pi*A)*np.amax(D)
    #print(RCS)
    RCS_db = 10*np.log10(RCS)
    print("encoding: " + str(opt_code))
    #plt.imshow(abs(D))
    print("RCS: " + str(RCS_db))    #print("D denominator: " + str(D_den[0]))
    return RCS_db
    #RCS_db_fit = 10*np.log10(32.08*pow(wavelength,2.187)/(4*pi*pow(A,1.0935)))
    #print(RCS_db_fit)



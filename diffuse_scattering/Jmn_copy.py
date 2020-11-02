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
    d=0.049
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

#def RCS_calc(opt_code,N):   
#    import numpy as np
#    from scipy import integrate
#    import matplotlib.pyplot as plt
#    from mpl_toolkits import mplot3d
#    from matplotlib.backends.backend_pdf import PdfPages
#    theta_i = np.pi/10
#    Eo = 1; 
#    k_x = 1;
#    k_y = 1;
#    d=0.049;
#    c=3e8;
#        #Number of elements in x-y directions in reflectarray
#    
#    wavelength = d
#    pi = np.pi
#    k = 2*pi/wavelength
#    
#    #N = 6
#    A=(N*d)**2
#    #opt_code = [0,0,1,0,1,1]  
#    
#    def F(t,p):
#        arr_factor = 0
#        for m in range(N):
#            for n in range(N):
#                arr_factor +=  gamma[m,n]*np.exp(-1j*d*(m*np.sin(t)*np.cos(p)+n*np.sin(t)*np.sin(p)))
#        return arr_factor
#
#    def den(t,p):
#        arr_factor = 0
#        for m in range(N):
#            for n in range(N):
#                arr_factor+=  gamma[m,n]*np.exp(-1j*d*(m*np.sin(t)*np.cos(p)+n*np.sin(t)*np.sin(p)))
#        return abs(arr_factor)**2*np.sin(t)*Es_theta_squared(t,p)
#    
#    Es_angular_y = np.zeros(shape=(100,100),dtype=np.complex_)
#    Es_angular_z = np.zeros(shape=(100,100),dtype=np.complex_)    
#    D = np.zeros(shape=(100,100),dtype=np.complex_)   
#    #F = np.zeros(shape=(100,100),dtype=np.complex_) 
#    E_freq = np.zeros(shape=(N,N),dtype=np.complex_)
#    phase = np.zeros((N,N))
#    
#    #matrix = np.repeat([opt_code],N,axis=0)
#   # for i in range(N):
#    #        if(matrix[0,i] == 1):
#   #             matrix[i] = ~matrix[i] +2
#    matrix=opt_code
#    for i in range(N):
#        for j in range(N):
#            phase[i,j] = matrix[i,j]
#    gamma = 2*phase-1
#    
#    #plt.imshow(phase)
#    #plt.show()
#    Es_freq = np.zeros((N,N))
#    #for x in range(N):
#    #    for y in range(N):
#    #        E_temp = np.zeros(shape=(N,N),dtype=np.complex_)
#    #        E_temp[x,y] = phase[10*x:10*x+10,10*y:10*y+10]*gamma[10*x:10*x+10,10*y:10*y+10]
#    #        Es_freq = Es_freq + np.fft.fftshift((np.fft.fft2(E_temp)))
#    
#    #plt.imshow(abs((Es_freq)))
#    #plt.show()
#    
#    Esy = lambda t,p: Eo*d**2*np.sinc(k*np.sin(t)*np.cos(p)*d/2)*np.sinc((k*np.sin(t)*np.sin(p)-np.sin(theta_i))*d/2)
#    Es_theta_y = lambda t, p: Esy(t,p)*(np.cos(t))
#    Es_theta_z = lambda t,p: Esy(t,p)*(np.sin(t)*np.sin(p))
#    Es_theta_squared = lambda t,p: abs(np.sqrt(Es_theta_y(t,p)**2+Es_theta_z(t,p)**2))**2
#    D_denominator = lambda t,p : np.sin(t)*Es_theta_squared(t,p)*abs(array_factor(k*np.sin(t)*np.cos(p),k*np.sin(t)*np.sin(p),gamma))**2
#    D_numerator = lambda t,p: 4*pi*Es_theta_squared(t,p)*abs(array_factor(t,p))**2
#    
#    D_den = integrate.dblquad(den,0,2*pi,0,pi/2,epsabs=1e-5)
#    arr_factor = np.zeros(shape=(100,100),dtype=np.complex_) 
#    f_y = np.zeros(shape=(100,100),dtype=np.complex_) 
#    f_z = np.zeros(shape=(100,100),dtype=np.complex_) 
#    Es = np.zeros(shape=(100,100),dtype=np.complex_) 
#    theta = np.linspace(0,pi/2,100)
#    phi = np.linspace(0,2*pi,100)
#    i=0
#    j=0
#    for t in theta:
#        for p in phi:
#            arr_factor[i,j]=F(t,p) 
#            Es[i,j] = Esy(t,p)
#            f_y[i,j] = Es_theta_y(t,p)
#            f_z[i,j] = Es_theta_z(t,p)
#            D[i,j] = (f_y[i,j]**2+f_z[i,j]**2)*abs(arr_factor[i,j])**2/D_den[0]
#            i+=1
#        j+=1
#        i=0
#    
#    RCS = wavelength**2/(4*pi*A)*np.amax(D)
#    #print(RCS)
#    RCS_db = 10*np.log10(RCS)
#    print("encoding: " + str(opt_code))
#    #plt.imshow(abs(D))
#    print("RCS: " + str(RCS_db))
#    print("standard deviation" + str(np.std(D)))
#    return RCS_db


import time
import matplotlib.pyplot as plt
N = 3
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
start_time = time.time()
sampleset = dimod.ExactSolver().sample_ising(h, J)
print(str(N) + " --- %s seconds ---" % (time.time() - start_time))

print(sampleset.first)
mat_dict = sampleset.first[0]
opt_code_2d = np.zeros(shape=(N,N))
for m in range(N):
    for n in range(N):
        dict_string = 'x'+str(m*N+n)
        opt_code_2d[m,n] = mat_dict[dict_string]
plt.imshow(opt_code_2d)        
print(J)

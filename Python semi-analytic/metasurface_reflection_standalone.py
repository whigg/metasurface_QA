def array_factor(t,p,gamma):
    arr_factor = 0
    for x in range(N):
        for y in range(N):
            arr_factor+= gamma[10*x,10*y] * np.exp(-1j*d*(x*np.sin(t)*np.cos(p)+y*np.sin(t)*np.sin(p)))
    return arr_factor


Eo = 1; 
k_x = 1;
k_y = 1;
d=0.049;
c=3e8;
    #Number of elements in x-y directions in reflectarray
 
wavelength = 0.049
pi = np.pi
k = 2*pi/wavelength

N = 7 
A=(N*d)**2
opt_code = [0,0,1,1,0,1,0]  

Es_angular_y = np.zeros(shape=(10*N,10*N),dtype=np.complex_)
Es_angular_z = np.zeros(shape=(10*N,10*N),dtype=np.complex_)    
D = np.zeros(shape=(10*N,10*N),dtype=np.complex_)   
f = np.zeros(shape=(10*N,10*N),dtype=np.complex_) 
F = np.zeros(shape=(10*N,10*N),dtype=np.complex_) 
E_freq = np.zeros(shape=(10*N,10*N),dtype=np.complex_)
phase = np.zeros((N*10,N*10))

matrix = np.repeat([opt_code],N,axis=0)
for i in range(N):
    if(matrix[0,i] == 1):
        matrix[i] = ~matrix[i] +2

for i in range(N):
    for j in range(N):
        phase[10*i:10*i+10,10*j:10*j+10] = matrix[i,j]
gamma = 2*phase-1

#plt.imshow(phase)
#plt.show()
Es_freq = np.zeros((10*N,10*N))
for x in range(N):
    for y in range(N):
        E_temp = np.zeros(shape=(10*N,10*N),dtype=np.complex_)
        E_temp[10*x:10*x+10,10*y:10*y+10] = phase[10*x:10*x+10,10*y:10*y+10]*gamma[10*x:10*x+10,10*y:10*y+10]
        Es_freq = Es_freq + np.fft.fftshift((np.fft.fft2(E_temp)))

#plt.imshow(abs((Es_freq)))
#plt.show()

Esy = lambda t,p: Eo*d**2*np.sinc(k*np.sin(t)*np.cos(p)*d/2)*np.sinc(k*np.sin(t)*np.sin(p)*d/2)
Es_theta_y = lambda t, p: Esy(t,p)*(np.cos(t))
Es_theta_z = lambda t,p: Esy(t,p)*(np.sin(t)*np.sin(p))
Es_theta_squared = lambda t,p: abs(np.sqrt(Es_theta_y(t,p)**2+Es_theta_z(t,p)**2))**2
D_denominator = lambda t,p : np.sin(t)*Es_theta_squared(t,p)*abs(array_factor(k*np.sin(t)*np.cos(p),k*np.sin(t)*np.sin(p),gamma))**2
D_numerator = lambda t,p: 4*pi*Es_theta_squared(t,p)*abs(array_factor(k*np.sin(t)*np.cos(p),k*np.sin(t)*np.sin(p),gamma))**2

D_den = integrate.dblquad(D_denominator,0,2*pi,0,pi/2,epsabs=1e-4)
arr_factor = np.complex(0)
theta = np.linspace(0,pi/2,N*10)
phi = np.linspace(0,2*pi,N*10)
i=0
j=0
for t in theta:
    for p in phi:
        f[i,j] = Es_theta_squared(t,p)
        Es_angular_y[i,j] = Es_theta_y(t,p)
        Es_angular_z[i,j] = Es_theta_z(t,p)
        F[i,j] = array_factor(t,p,gamma)
        E_freq[i,j] = Esy(t,p)
        D[i,j] = D_numerator(t,p)/D_den[0]
        i+=1
    j+=1
    i=0
    
#plt.imshow(D.real)
RCS = wavelength**2/(4*N**2*d)*np.amax(D)
#print(RCS)
RCS_db = 10*np.log10(RCS)
#print("encoding: " + str(opt_code))
print("RCS: " + str(RCS_db))
#print("D denominator: " + str(D_den[0]))

#RCS_db_fit = 10*np.log10(32.08*pow(wavelength,2.187)/(4*pi*pow(A,1.0935)))
#print(RCS_db_fit)

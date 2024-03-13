import json

import numpy as np
import matplotlib.pyplot as plt

from scipy import stats
from scipy.optimize import curve_fit

kB = 1.0

with open('./results/simulation.json') as f:
    sim = json.load(f)

T = sim['global']['ensemble']['data'][0][sim['global']['ensemble']['labels'].index('temperature')]
print('Temperature: {}'.format(T))

name,mass,radius,charge = sim['global']['types']['data'][0]
print('Particle type, name: {}, mass: {}, radius: {}, charge: {}'.format(name,mass,radius,charge))

N = len(sim['state']['data'])

print('Number of particles: {}'.format(N))

####################################

print("Testing EulerMaruyama...")

dtEulerMaruyama = sim['integrator']['eulerMaruyama']['parameters']['timeStep']

#Check EULER MARUYAMA, MSD is compareted with the analytical solution
v = sim['integrator']['eulerMaruyama']['parameters']['viscosity']

msdEulerMaruyama = np.loadtxt('./results/msdEulerMaruyama.dat',skiprows=1)
time      = (msdEulerMaruyama[:,0]-msdEulerMaruyama[0,0])*dtEulerMaruyama
msdSim    = msdEulerMaruyama[:,1]

#Fit the MSD with a linear function
slope, intercept, r_value, p_value, std_err = stats.linregress(time,msdSim)
#msd = 2*D*d*t # Where D is the diffusion coefficient, d is the dimensionality and t is the time
Dsim  = slope/(2.0*3.0)

M = 1.0/(6.0*np.pi*v*radius)
Dana = kB*T*M
print('D trans EulerMaruyama: {} (sim) vs {} (ana)'.format(Dsim,Dana))

#Plot the MSD
plt.plot(time,msdSim,'o',label='Simulation')
plt.plot(time,2.0*Dana*3.0*time,'-',label='Analytical')
plt.xlabel('Time')
plt.ylabel('MSD')
plt.title('EulerMaruyama test (translational diffusion)')
plt.legend()
plt.show()

####################################

print("Testing EulerMaruyamaRigid...")

dtEulerMaruyamaRigid = sim['integrator']['eulerMaruyamaRigid']['parameters']['timeStep']

#Check EULER MARUYAMA, MSD is compareted with the analytical solution
v = sim['integrator']['eulerMaruyamaRigid']['parameters']['viscosity']

msdEulerMaruyamaRigid = np.loadtxt('./results/msdEulerMaruyamaRigid.dat',skiprows=1)
time      = (msdEulerMaruyamaRigid[:,0]-msdEulerMaruyamaRigid[0,0])*dtEulerMaruyamaRigid
msdSim    = msdEulerMaruyamaRigid[:,1]

#Fit the MSD with a linear function
slope, intercept, r_value, p_value, std_err = stats.linregress(time,msdSim)
#msd = 2*D*d*t # Where D is the diffusion coefficient, d is the dimensionality and t is the time
Dsim  = slope/(2.0*3.0)

M = 1.0/(6.0*np.pi*v*radius)
Dana = kB*T*M
print('D trans EulerMaruyamaRigid: {} (sim) vs {} (ana)'.format(Dsim,Dana))

#Plot the MSD
plt.plot(time,msdSim,'o',label='Simulation')
plt.plot(time,2.0*Dana*3.0*time,'-',label='Analytical')
plt.xlabel('Time')
plt.ylabel('MSD')
plt.title('EulerMaruyamaRigid test (translational diffusion)')
plt.legend()
plt.show()

#Checkin rotational diffusion

macEulerMaruyamaRigid = np.loadtxt('./results/macEulerMaruyamaRigid.dat',skiprows=1)
time     = (macEulerMaruyamaRigid[:,0]-macEulerMaruyamaRigid[0,0])*dtEulerMaruyamaRigid
macSim   = macEulerMaruyamaRigid[:,1]

#Fit mean angular correlation (MAC) with a exponential decay

def func(t, D):
    return np.exp(-2.0*D*t)

popt, pcov = curve_fit(func, time, macSim)
Dsim = popt[0]

M = 1.0/(8.0*np.pi*v*radius**3)
Dana = kB*T*M
print('D rot EulerMaruyamaRigid: {} (sim) vs {} (ana)'.format(Dsim,Dana))

#Plot the MAC
plt.plot(time,macSim,'o',label='Simulation')
plt.plot(time,func(time,Dana),'-',label='Analytical')
plt.xlabel('Time')
plt.ylabel('MAC')
plt.title('EulerMaruyamaRigid test (rotational diffusion)')
plt.legend()
plt.show()

####################################

def maxwellVelocityDistribution(v,m,kBT):
    return np.power(m/(2.0*np.pi*kBT),1.5)*4.0*np.pi*v*v*np.exp(-m*v*v/(2.0*kBT))

####################################

print("Testing BBK...")
dtBBK = sim['integrator']['bbk']['parameters']['timeStep']

#Check BBK velocity distribution

vBBK = np.loadtxt('./results/BBK.vel')
vBBK = vBBK[:,3] #The last column stores the velocity module

#Remove the first N lines
vBBK = vBBK[N:]

#Plot histogram and compare with the Maxwell-Boltzmann distribution
plt.hist(vBBK,bins=100,density=True,label='Simulation')
v = np.linspace(np.min(vBBK),np.max(vBBK),1000)
plt.plot(v,maxwellVelocityDistribution(v,mass,kB*T),label='Maxwell-Boltzmann')
plt.xlabel('Velocity')
plt.ylabel('Probability')
plt.title('BBK test (velocity distribution)')
plt.legend()
plt.show()

print("Testing GJF...")

dtGJF = sim['integrator']['gjf']['parameters']['timeStep']

#Check GJF velocity distribution

vGJF = np.loadtxt('./results/GJF.vel')
vGJF = vGJF[:,3] #The last column stores the velocity module

#Remove the first N lines
vGJF = vGJF[N:]

#Plot histogram and compare with the Maxwell-Boltzmann distribution
plt.hist(vGJF,bins=100,density=True,label='Simulation')
v = np.linspace(np.min(vGJF),np.max(vGJF),1000)
plt.plot(v,maxwellVelocityDistribution(v,mass,kB*T),label='Maxwell-Boltzmann')
plt.xlabel('Velocity')
plt.ylabel('Probability')
plt.title('GJF test (velocity distribution)')
plt.legend()
plt.show()






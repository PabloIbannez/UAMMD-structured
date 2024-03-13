from scipy.integrate import odeint
from numpy import pi, zeros, linspace, sin
import numpy as np
from scipy.signal import find_peaks

#----------------------------------------utils----------------------------------------------------#
def getLastCycle(field, magnetization, pointsPerCycle=-1, ncycles=-1):
     if (pointsPerCycle>0 and ncycles>0):
          if (pointsPerCycle*ncycles + 1 == len(field)):
               return field[-pointsPerCycle-1:], magnetization[-pointsPerCycle-1:]
          else:
               return field[pointsPerCycle*ncycles-1:], magnetization[pointsPerCycle*ncycles-1:]

     else:
          lastField = field[-1]
          distanceToLastField = abs(field[:-2]-lastField)
          minimumDistances, _ = find_peaks(-distanceToLastField)
          return field[minimumDistances[-2]:], magnetization[minimumDistances[-2]:]


     
#-------------------------------Solvers FP equation-----------------------------------------------#

def compute_dadt(a,t,tneel,alpha0,alphak,w):
    nindex = len(a)
    field = alpha0*sin(t*w)
    dan = zeros(nindex)
    dan[0] = 0
    for n in range(1,nindex-2):
        termField = (-a[n]+field*(a[n-1]/(2*n-1)-a[n+1]/(2*n+3)))
        termAnisotropy = alphak*(n*a[n]/((2*n-1)*(2*n+1))-(n+1)*a[n]/((2*n+1)*(2*n+3))-(n+2)*a[n+2]/((2*n+3)*(2*n+5)))
        if n>1:
            termAnisotropy += alphak*(n-1)*a[n-2]/((2*n-3)*(2*n-1))
        dan[n]=0.5/tneel*n*(n+1)*(termField+termAnisotropy)
    return dan


def integrateFP(tn, alpha0, alphak, w, nindex, time):
    an_0 = zeros(nindex)
    an_0[0] = 0.5
    an_t = odeint(compute_dadt, an_0, time, args = (tn,alpha0, alphak, w))
    return an_t


def computeMagnetization(kbT, rc, K, damp, m0, gamma0, f, b0,
                              nindex = 30, ncycles = 10, pointsPerCycle = 500, lastCycle = None):

    t = m0*(1+damp*damp)/(2*kbT*gamma0*damp)
    alpha0 = m0*b0/kbT
    Vc = 4./3.*np.pi*rc*rc*rc
    alphak = 2*K*Vc/kbT
    numberPoints = int(ncycles*pointsPerCycle)+1
    time = linspace(0,ncycles/f, numberPoints)
    if (lastCycle is not None):
         lastCycle+=time[-1]
         time = np.concatenate((time,lastCycle))
    an_t = integrateFP(t, alpha0, alphak, 2*pi*f, nindex, time)
    magnetization = 2*an_t[:,1]/3.0
    field = b0*sin(2*pi*f*time)
    return time, field, magnetization


def computeCycle(kbT, rc, K, damp, m0, gamma0, f, b0, nindex = 30,
                      ncycles = 10, pointsPerCycle = 500, lastCycle = None):
     time, field, magnetization = computeMagnetization(kbT, rc, K, damp, m0, gamma0, f, b0,
                                                            nindex, ncycles, pointsPerCycle,
                                                            lastCycle)
     field, magnetization = getLastCycle(field, magnetization, pointsPerCycle, ncycles)
     return field, magnetization       

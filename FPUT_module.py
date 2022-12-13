import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import animation
from matplotlib.animation import FuncAnimation, ArtistAnimation

import scipy as sp
import scipy.constants as spc

from IPython.core.display import HTML
from IPython.display import display, Math, Latex

import os, time

def customPlt():
    plt.rcParams.update({
                  'text.usetex':True,
                  'font.family':'serif',
              'font.sans-serif':'Times New Roman',
             'mathtext.fontset':'cm',
               'axes.edgecolor':'gray',
               'axes.axisbelow':True,
               'axes.labelsize':18,
               'axes.titlesize':20,
        'legend.title_fontsize':18,
              'legend.fontsize':16,
                    'font.size':18,
             'figure.titlesize':22,
              'xtick.labelsize':16,
              'ytick.labelsize':16,
              'xtick.direction':'out',
              'ytick.direction':'out',
                  'xtick.color':'black',
                  'ytick.color':'black',
                  'axes3d.grid':False,
        'animation.embed_limit':2**128
    })

pi = spc.pi
amu_kg = spc.atomic_mass # [kg]
eV = spc.e # [J]
hbar = spc.hbar # [J/s]


def fixed_fput(q0, p0, N = 32, epsilon = 0, dt = 0.1, steps = 100, method = 'leapfrog'):
    """ Solve for the dynamics of the discrete alpha-FPUT lattice.
        Args:
            - q0(np.array) : initial displacement profile.
            - p0(np.array) : initial momentum profile.
        Kwargs:
            - N (int) : Number of active particles in the lattice (excluding the fixed ends).
            - epsilon (float) : non-linear coupling strength of cubic potential.
            - dt (float) : time step.
            - steps (int) : Number of time steps.
            - method (str) : the numerical integration method used to solve for the dynamical variables q and p.
                - User can choose any one of the following methods: 
                    1. Forward Euler ('fwd_Eul').
                    2. Symplectic Euler ('sym_Eul').
                    3. Leapfrog ('leapfrog').
        Returns:
            t, q, p (dimensionless arrays for simulation time, displacement profiles and momentum profiles).
    """
    def F(x):
        d1 = np.roll(x,-1)+np.roll(x,1)-2*x
        d2 = np.roll(x,-1)-np.roll(x,1)
        d1[0] = d1[-1] = d2[0] = d2[-1] = 0.0        
        force = d1*(1 + 0.5*epsilon*d2)
        force[0] = force[-1] = 0.0
        return force
    
    t = dt*np.arange(steps)
    q = np.zeros((steps+2, N+2))
    p = np.zeros((steps+2, N+2))
    
    q[0] = q[1] = q0
    p[0] = p[1] = p0

    for i in range(1,steps+1):
        q[:,0] = q[:,-1] = p[:,0] = p[:,-1] = 0.0        
        if (method == 'fwd_Eul'):
            p[i+1] = p[i] + dt*F(q[i])
            q[i+1] = q[i] + dt*p[i]            
        if (method == 'sym_Eul'):
            p[i+1] = p[i] + dt*F(q[i])
            q[i+1] = q[i] + dt*p[i+1]
        if (method == 'leapfrog'):
            p_mid = p[i] + 0.5*dt*F(q[i])
            q[i+1] = q[i] + dt*p_mid
            p[i+1] = p_mid + dt*F(q[i+1])
        q[:,0] = q[:,-1] = p[:,0] = p[:,-1] = 0.0
    t = t    
    q = q[1:-1]
    p = p[1:-1]
    return t, q, p

def initMode(lattice, mode):
    N = len(lattice)
    A = np.sqrt(2/(N-1))
    return A * np.sin(pi*mode*lattice/(N-1)) 

def dst(q):
    N = q.shape[1]
    A = np.sqrt(2/(N-1))
    k = np.arange(N)
    n = k.reshape((N,1))
    sin_k = np.sin(pi*k*n/(N-1))
    Q = A*np.dot(q, sin_k)
    return Q

def getEnergy(Q, P):
    Qmode,Pmode = Q.T, P.T
    N = Q.shape[1]
    E_k = np.zeros_like(Qmode)
    for k in range(N):
        freq_k = 2*np.sin(0.5*k*pi/(N-1))
        E_k[k] = (Pmode[k]**2 + freq_k**2 * Qmode[k]**2)/2
    return E_k.T

def kdv_FTCS(f_next, f_n, dz, dt, delta, order):
    
    c1 = dt/dz
    c2 = (delta/dz)**2
    Nz = len(f_next)
    f_next[0] = f_next[Nz-1] = f_n[0] = f_n[Nz-1]
    if order == 2:
        for m in range(2, Nz-2):
            d1 = f_n[m+1]-f_n[m-1]
            d2 = f_n[m+2] - f_n[m-2]            
            f_next[m] = f_n[m] - c1*(f_n[m]*(0.5*d1) + c2*(0.5*d2-d1))
    elif order == 6:
        for m in range(4, Nz-4):
            d1 = f_n[m+1] - f_n[m-1]
            d2 = f_n[m+2] - f_n[m-2]
            d3 = f_n[m+3] - f_n[m-3]
            d4 = f_n[m+4] - f_n[m-4]
            f_next[m] = f_n[m] - c1*(f_n[m]*(d3*1/60 - d2*3/20 + d1*3/4) + c2*(d4*7/240 - d3*3/10 + d2*169/120 - d1*61/30))    
    f_next[0] = f_next[Nz-1] = f_n[0] = f_n[Nz-1]
    
def kdv_solveFTCS(f0, z, dz, t, dt, delta, order = 6):
    f_next = np.zeros(len(z))
    f = []
    f.append(f0.copy())
    for n in range(len(t)-1):
        kdv_FTCS(f_next, f[n], dz, dt, delta, order)
        f.append(f_next.copy())
    return np.array(f)

def kdv_LAX(f_next, f_n, dz, dt, delta):
    c1 = dt/dz
    c2 = (delta/dz)**2
    Nz = len(f_next)
    f_next[0] = f_next[Nz-1] = f_n[0] = f_n[Nz-1]
    for m in range(4, Nz-4):
        d1 = f_n[m+1] - f_n[m-1]
        d2 = f_n[m+2] - f_n[m-2]
        d3 = f_n[m+3] - f_n[m-3]
        d4 = f_n[m+4] - f_n[m-4]
        lax1 = (1/2)*(f_n[m+1] - f_n[m-1])
        lax2 = (1/4)*((f_n[m+1])**2 - (f_n[m-1])**2)
        f_next[m] = f_n[m] + lax1 - c1*(lax2 + c2*(d4*7/240 - d3*3/10 + d2*169/120 - d1*61/30))  
    f_next[0] = f_next[Nz-1] = f_n[0] = f_n[Nz-1]    
    
def kdv_solveLAX(f0, z, dz, t, dt, delta):
    f_next = np.zeros(len(z))
    f = []
    f.append(f0.copy())
    for n in range(len(t)-1):
        kdv_LAX(f_next, f[n], dz, dt, delta)
        f.append(f_next.copy())
    return np.array(f)

def gaussIC(z, z0, sigma):
    A = np.sqrt(0.5/np.pi)/sigma
    p = -0.5 * ((z - z0)/sigma)**2
    return A*np.exp(p)

def solitonIC(z,z0,v):
    f = np.sqrt(v/2) * np.cosh( np.sqrt(v/2)*(z - z0) )
    return f**(-2)

def check_stability(dt, dz, delta):
    G = dt/(dz**3)
    Gmax = 1.5*np.sqrt(3)/(delta**2)
    if G <= Gmax:
        display(Latex(r'$\mathrm{Stability\,Condition\,Met}$!'))
        display(Latex(r'$\displaystyle \frac{\Delta \tau}{\Delta z^3} = $ ' 
                      + r'${:f} \leq $'.format(G) + r'$\displaystyle\frac{3\sqrt{3}}{2\delta^2}$'))
        return True
    else:
        display(Latex(r'Heavily unstable! Try a smaller time step?'))
        return False

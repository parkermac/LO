"""
Code to test whether I can recover the incident and reflected waves
from a case where I create a synthetic signal with known incident and
reflected waves and known friction.
"""

import numpy as np
import cmath
import utide
from datetime import datetime, timedelta
import matplotlib.pyplot as plt

om = 2*np.pi/12.420601202671868 # frequency, radians per hour: M2
th = np.arange(0,24*365)    # time in hours

def get_hm(th, E):
    hm = utide.solve(th/24, E, lat=0, constit=['M2'], nodal=False, trend=False, phase='raw')
    A = hm.A[0]
    G = np.pi*hm.g[0]/180
    return A, G

# zero phase shift signal to get phase reference gr0
eta0 = np.cos(-om*th)
A0, G0 = get_hm(th, eta0)

if False:
    # A test of the methods used here. RESULT: it works perfectly
    # synthetic signal to analyze
    eta = np.cos(-om*th - np.pi/8)
    A, G = get_hm(th, eta)
    # reconstruct the signal
    eta_r = A * np.cos(-om*th + G - G0)
    # PLOTTING
    td = th/24
    plt.close('all')
    fig = plt.figure(figsize=(16,8))
    ax = fig.add_subplot(111)
    ax.plot(td, eta, '-c', lw=5)
    ax.plot(td, eta_r, '--b', lw=2)
    ax.set_xlim(0,3)
    ax.grid(True)
    ax.set_xlabel('Time (days)')
    plt.show()
    
# create incident (+, p) and reflected waves (-, m), with friction
g = 9.8
H=40
fric = 1 # R/omega
alpha = np.sqrt(g/H)/np.sqrt(1 + fric*1j)

# Incident
ep = 1 * np.cos(-om*th)
Aep, Gep = get_hm(th, ep)
Cep = cmath.rect(Aep, Gep)  # complex amplitude of eta
Cvp = alpha * Cep           # complex amplitude of v

if False:
    # Check on complex reconstruction. RESULT: it works perfectly
    ep_r0 = Aep * np.cos(-om*th + Gep - G0)
    ep_r1 = (Cep * np.exp(-th*om*1j - G0*1j)).real
    # PLOTTING
    td = th/24
    plt.close('all')
    fig = plt.figure(figsize=(16,8))
    ax = fig.add_subplot(111)
    ax.plot(td, ep_r0, '-c', lw=5)
    ax.plot(td, ep_r1, '--b', lw=2)
    ax.plot(td, ep, '*', c='purple', ms=10)
    ax.set_xlim(0,3)
    ax.grid(True)
    ax.set_xlabel('Time (days)')
    plt.show()

# Reflected
em = .5 * np.cos(-om*th + np.pi/8)
Aem, Gem = get_hm(th, em)
Cem = cmath.rect(Aem, Gem)  # complex amplitude of eta
Cvm = alpha * Cem           # complex amplitude of v

# combined signal
Ce = Cep + Cem
Cv = Cvp - Cvm
ee = (Ce * np.exp(-th*om*1j - G0*1j)).real
vv = (Cv * np.exp(-th*om*1j - G0*1j)).real
Aee, Gee = get_hm(th, ee)
Avv, Gvv = get_hm(th, vv)
Ce_alt = cmath.rect(Aee, Gee)  # complex amplitude of eta
Cv_alt = cmath.rect(Avv, Gvv)  # complex amplitude of v

# recreate the incident and reflected waves: works perfectly
Cep_alt = (Ce_alt + Cv_alt/alpha)/2
Cem_alt = (Ce_alt - Cv_alt/alpha)/2

# what if I don't know alpha but instead guess?


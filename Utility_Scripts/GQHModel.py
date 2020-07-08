# GQ/H function
from math import sqrt
import numpy as np

# Genralized Quadratic Hyperboloic Model
def GQH(tmax, G0, t1, t2, t3, t4, t5, gamma):
    gr = tmax / G0
    tau = []
    for g in gamma:
        if g == 0.0:
            tau.append(0.0)
            continue
        tt = t1 + t2 * (t4*(g/gr)**t5) / (t3**t5 + t4*(g/gr)**t5)
        try:
            thisTau = 0.5*tmax / tt * (1 + g/gr - sqrt(((1+(g/gr))**2 - 4*tt*g/gr)))
        except:
            thisTau = 0.0
        tau.append(thisTau)
    return np.array(tau)

# Calculate hysteretic loops based on Masing rule
def GQH_Masing(tmax, G0, t1, t2, t3, t4, t5, gamma):
    loops_gamma = []
    loops_tau   = []
    for ii in range(len(gamma)):
        g = np.linspace(0, gamma[ii], 200)
        t = GQH(tmax, G0, t1, t2, t3, t4, t5, g)
        g = 2.0 * g - g[-1]
        t = 2.0 * t - t[-1]
        g = np.append(g, np.flip(g))
        t = np.append(t, -t)

        loops_gamma.append(g)
        loops_tau.append(t)
    return loops_gamma, loops_tau

# calculate the damping curve for GQH model
def GQH_damping(tmax, G0, t1, t2, t3, t4, t5, gamma):
    loops_gamma, loops_tau = GQH_Masing(tmax, G0, t1, t2, t3, t4, t5, gamma)
    damp = []
    for ii in range(len(loops_gamma)):
        d_gamma = np.diff(loops_gamma[ii])
        tau_med = loops_tau[ii][0:-1] + 0.5 * np.diff(loops_tau[ii])
        loop_area = np.dot(d_gamma, tau_med)
        maxGam_index = np.argmax(loops_gamma[ii])
        elastic_energy = 0.5 * loops_gamma[ii][maxGam_index] * loops_tau[ii][maxGam_index]
        damp.append(0.25 * loop_area / elastic_energy / np.pi)
    return np.array(damp)
    
# return the G/Gmax curve for the GQH model
def GQH_modulusreduction(tmax, G0, t1, t2, t3, t4, t5, gamma):
    tau = GQH(tmax, G0, t1, t2, t3, t4, t5, gamma)
    return tau/gamma/G0
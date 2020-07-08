# calibrate GQH
import numpy as np
import matplotlib.pyplot as plt

from GQHModel import GQH_damping, GQH_modulusreduction
from ModulusReduction import Darendeli, Menq, VuceticDobry, Zhang2005
from scipy.optimize import minimize, least_squares, fmin_slsqp

# objective function for minimization (G/Gmax)
def errorGGmax(matProps, tau_max, G_max, gamma, GGmax):
    # tmax, G0, t1, t2, t3, t4, t5
    GQH_GGmax = GQH_modulusreduction(tau_max, G_max, matProps[0], matProps[1], matProps[2], 1.0, matProps[3], gamma)
    return np.dot(GGmax[np.where(gamma <= 0.0005)] - GQH_GGmax[np.where(gamma <= 0.0005)], GGmax[np.where(gamma <= 0.0005)] - GQH_GGmax[np.where(gamma <= 0.0005)])

# objective function for minimization (Damping)
def errorDamp(factProps, t1, t2, t3, t4, t5, tau_max, G_max, gamma, damp):
    GQH_damp  = GQH_damping(tau_max, G_max, t1, t2, t3, t4, t5, gamma)
    GQH_GGmax = GQH_modulusreduction(tau_max, G_max, t1, t2, t3, t4, t5, gamma)
    factor = factProps[0] - factProps[1] * (1.0 - GQH_GGmax)**factProps[2]
    GQH_damp  = factProps[3] + factor * GQH_damp
    return np.dot((damp -GQH_damp)/damp, (damp - GQH_damp)/damp)

# residual function for the GGmax least squares method
def GHQ_tt_res(t, tt_GGmax, tau_max, G_max, gamma):
    gr = tau_max/G_max
    tt = t[0]+t[1]*((gamma/gr)**t[3])/(t[2]**t[3]+(gamma/gr)**t[3])
    res = tt - tt_GGmax
    return res

# residual function for the GGmax least squares when strength is not captured (step 1)
def GHQ_tt_res95_1(t, tt_GGmax, tau_max, G_max, gamma, t3, t5):
    gamma = np.append(gamma, [0.1])
    gr = tau_max/G_max
    tt95 = (0.95 + 0.95 * 0.1 / gr - 0.1 / gr) / (0.95**2.0)
    t1 = t[0]
    t2 = t[1]
    tt = t1+t2*((gamma/gr)**t5)/(t3**t5+(gamma/gr)**t5)
    res = tt[0:-1] - tt_GGmax
    res = np.append(res, [100.0*(tt[-1] - tt95)])
    return res

# residual function for the GGmax least squares when strength is not captured (step 2)
def GHQ_tt_res95_2(t, tt_GGmax, tau_max, G_max, gamma, t1, t2):
    gamma = np.append(gamma, [0.1])
    gr = tau_max/G_max
    tt95 = (0.95 + 0.95 * 0.1 / gr - 0.1 / gr) / (0.95**2.0)
    t3 = t[0]
    t5 = t[1]
    tt = t1+t2*((gamma/gr)**t5)/(t3**t5+(gamma/gr)**t5)
    res = tt[0:-1] - tt_GGmax
    res = np.append(res, [100.0*(tt[-1] - tt95)])
    return res

# residual function for the damping least squares method
def resDamp(factProps, t1, t2, t3, t4, t5, tau_max, G_max, gamma, damp):
    GQH_damp  = GQH_damping(tau_max, G_max, t1, t2, t3, t4, t5, gamma)
    GQH_GGmax = GQH_modulusreduction(tau_max, G_max, t1, t2, t3, t4, t5, gamma)
    GQH_GGmax[np.where(GQH_GGmax > 1.0)] = 1.0
    factor = factProps[0] - factProps[1] * (1.0 - GQH_GGmax)**factProps[2]
    GQH_damp  = factor * GQH_damp
    return damp - GQH_damp

# constrain t1 + t2 < 1.0
def constraint1(t, tau_max, G_max, gamma, GGmax):
    return np.atleast_1d(1.0 - t[0] - t[1])

# constrain theta_tau < 1.0
def constraint2(t, tau_max, G_max, gamma, GGmax):
    gr = tau_max/G_max
    tt = t[0]+t[1]*((gamma/gr)**t[3])/(t[2]**t[3]+(gamma/gr)**t[3])
    return np.atleast_1d(np.sum(1 - tt))

# constrain theta_tau < 1.0
def constraint3(t, tau_max, G_max, gamma, GGmax):
    GQH_GGmax = GQH_modulusreduction(tau_max, G_max, t[0], t[1], t[2], 1.0, t[3], gamma)
    return np.atleast_1d(G_max*GQH_GGmax[-1]*gamma[-1] - 0.95*tau_max)

# calibration function for GQH model based on Darendeli curves
def calibrateGQH_Darendeli(p, PI, OCR, tau_max, G_max, frq=1, Ncyc=10, Patm=100.0):
    GGmax, damp, gamma = Darendeli(p, PI, OCR, frq, Ncyc, Patm)
    initProps = [0.0, 0.5, 0.0, 1.0]
    res = minimize(errorGGmax, initProps, (tau_max, G_max, gamma, GGmax), method="Nelder-Mead")
    t1 = res.x[0]
    t2 = res.x[1]
    t3 = res.x[2]
    t4 =  1.0
    t5 = res.x[3]
    
    initProps = [0.55,0.25,10.0, 0.0]
    res = minimize(errorDamp,initProps,(t1, t2, t3, t4, t5, tau_max, G_max, gamma, damp/100.0), method="Nelder-Mead")

    p1 = res.x[0]
    p2 = res.x[1]
    p3 = res.x[2]
    D  = res.x[3]

    return t1, t2, t3, t4, t5, p1, p2, p3, D

# calibration function for GQH model based on Menq curves
def calibrateGQH_Menq(p, Cu, D50, tau_max, G_max, Ncyc=10, Patm=100.0):
    GGmax, damp, gamma = Menq(p, Cu, D50, Ncyc, Patm)
    initProps = [0.0, 0.5, 0.0, 1.0]
    res = minimize(errorGGmax, initProps, (tau_max, G_max, gamma, GGmax), method="Nelder-Mead")
    t1 = res.x[0]
    t2 = res.x[1]
    t3 = res.x[2]
    t4 =  1.0
    t5 = res.x[3]
    
    initProps = [1.0,0.0,1.0, 0.0]
    res = minimize(errorDamp,initProps,(t1, t2, t3, t4, t5, tau_max, G_max, gamma, damp/100.0), method="Nelder-Mead")

    p1 = res.x[0]
    p2 = res.x[1]
    p3 = res.x[2]
    D  = res.x[3]

    return t1, t2, t3, t4, t5, p1, p2, p3, D

# calibration function for GQH model based on Darendeli curves - using least squares method
def calibrateGQH_Darendeli_LS(p, PI, OCR, tau_max, G_max, frq=1, Ncyc=10, Patm=100.0):
    # get the reference modulus reduction curve
    GGmax2, damp2, gamma2 = Darendeli(p, PI, OCR, frq, Ncyc, Patm)

    # set a cut-off strain for optimization
    eps_cutoff = 0.0005
    index = np.where(gamma2 > eps_cutoff)
    GGmax = np.delete(GGmax2, index)
    gamma = np.delete(gamma2, index)

    # calculate the theta_tau for the reference GGmax curve
    gr = tau_max/G_max
    Gam = (gamma/gr)
    tt = (GGmax + GGmax * Gam - 1.0) / (GGmax**2.0 * Gam)
    tt[np.where(tt > 1.0)] = 1.0

    # initial guess
    t10 = -1.0
    t20 = -3.0
    t30 = 1.0
    t50 = 0.5

    # perform the least squares analysis
    initProp = [t10,t20,t30,t50]
    # res = least_squares(GHQ_tt_res, initProp, args=(tt, tau_max, G_max, gamma), bounds=([-20.0, -20.0, 0.0001, 0.0],[20.0, 20.0, np.inf, 0.99]))
    # t1 = res.x[0]
    # t2 = res.x[1]
    # t3 = res.x[2]
    # t4 = 1.0
    # t5 = res.x[3]
    res = fmin_slsqp(errorGGmax, initProp, args=(tau_max, G_max, gamma2, GGmax2), bounds=[(-20,20),(-20,20),(0,np.inf),(0.1,0.99)], ieqcons=[constraint1, constraint2, constraint3], iprint=-1)
    t1 = res[0]
    t2 = res[1]
    t3 = res[2]
    t4 = 1.0
    t5 = res[3]

    # check the strength and the optimized values
    GQH_GGmax = GQH_modulusreduction(tau_max,G_max, t1, t2, t3, t4, t5, gamma2)
    # app_strength = (G_max*GQH_GGmax*gamma2)[-1]
    # if  (app_strength < 0.95 * tau_max) or (t1 + t2 > 1.0):
    #     initProp = [t10, t20]
    #     res = least_squares(GHQ_tt_res95_1, initProp, args=(tt, tau_max, G_max, gamma, t3, t5), bounds=([-20.0, -20.0],[20.0, 20.0]))
    #     t1 = res.x[0]
    #     t2 = res.x[1]
    #     initProp = [t30, t50]
    #     res = least_squares(GHQ_tt_res95_2, initProp, args=(tt, tau_max, G_max, gamma, t1, t2), bounds=([0.0001, 0.0],[np.inf, 0.99]))
    #     t3 = res.x[0]
    #     t5 = res.x[1]
    
    # optimize the damping parameters
    initProps = [1.0,0.0,1.0]
    res = least_squares(resDamp,initProps, args=(t1, t2, t3, t4, t5, tau_max, G_max, gamma2, damp2/100.0), bounds=([0.0, 0.0, 0.0],[np.inf, np.inf, np.inf]))
    p1 = res.x[0]
    p2 = res.x[1]
    p3 = res.x[2]

    # calculate the minimum damping and revise the parameters
    GQH_damp_org  = GQH_damping(tau_max, G_max, t1, t2, t3, t4, t5, gamma)
    factor = p1 - p2 * (1.0 - GQH_GGmax)**p3
    D = damp2[0]/100.0 - factor[0] * GQH_damp_org[0]
    res = least_squares(resDamp,initProps, args=(t1, t2, t3, t4, t5, tau_max, G_max, gamma2, damp2/100.0 - D), bounds=([0.0, 0.0, 0.0],[np.inf, np.inf, np.inf]))
    p1 = res.x[0]
    p2 = res.x[1]
    p3 = res.x[2]

    return t1, t2, t3, t4, t5, p1, p2, p3, D

# calibration function for GQH model based on Menq curves (least squares method)
def calibrateGQH_Menq_LS(p, Cu, D50, tau_max, G_max, Ncyc=10, Patm=100.0):
    GGmax2, damp2, gamma2 = Menq(p, Cu, D50, Ncyc, Patm)

    # set a cut-off strain for optimization
    eps_cutoff = 0.0005
    index = np.where(gamma2 > eps_cutoff)
    GGmax = np.delete(GGmax2, index)
    gamma = np.delete(gamma2, index)

    # calculate the theta_tau for the reference GGmax curve
    gr = tau_max/G_max
    Gam = (gamma/gr)
    tt = (GGmax + GGmax * Gam - 1.0) / (GGmax**2.0 * Gam)
    tt[np.where(tt > 1.0)] = 1.0

    # initial guess
    t10 = -1.0
    t20 = -3.0
    t30 = 1.0
    t50 = 0.5

    # perform the least squares analysis
    initProp = [t10,t20,t30,t50]
    # res = least_squares(GHQ_tt_res, initProp, args=(tt, tau_max, G_max, gamma), bounds=([-20.0, -20.0, 0.0001, 0.0],[20.0, 20.0, np.inf, 0.99]))
    # t1 = res.x[0]
    # t2 = res.x[1]
    # t3 = res.x[2]
    # t4 = 1.0
    # t5 = res.x[3]
    res = fmin_slsqp(errorGGmax, initProp, args=(tau_max, G_max, gamma2, GGmax2), bounds=[(-20,20),(-20,20),(0,np.inf),(0.1,0.99)], ieqcons=[constraint1, constraint2, constraint3], iprint=-1)
    t1 = res[0]
    t2 = res[1]
    t3 = res[2]
    t4 = 1.0
    t5 = res[3]

    # check the strength and the optimized values
    GQH_GGmax = GQH_modulusreduction(tau_max,G_max, t1, t2, t3, t4, t5, gamma2)
    # app_strength = (G_max*GQH_GGmax*gamma2)[-1]
    # if  (app_strength < 0.95 * tau_max) or (t1 + t2 > 1.0):
    #     initProp = [t10, t20]
    #     res = least_squares(GHQ_tt_res95_1, initProp, args=(tt, tau_max, G_max, gamma, t3, t5), bounds=([-20.0, -20.0],[20.0, 20.0]))
    #     t1 = res.x[0]
    #     t2 = res.x[1]
    #     initProp = [t30, t50]
    #     res = least_squares(GHQ_tt_res95_2, initProp, args=(tt, tau_max, G_max, gamma, t1, t2), bounds=([0.0001, 0.0],[np.inf, 0.99]))
    #     t3 = res.x[0]
    #     t5 = res.x[1]
    
    # optimize the damping parameters
    initProps = [1.0,0.0,1.0]
    res = least_squares(resDamp,initProps, args=(t1, t2, t3, t4, t5, tau_max, G_max, gamma2, damp2/100.0), bounds=([0.0, 0.0, 0.0],[np.inf, np.inf, np.inf]))
    p1 = res.x[0]
    p2 = res.x[1]
    p3 = res.x[2]

    # calculate the minimum damping and revise the parameters
    GQH_damp_org  = GQH_damping(tau_max, G_max, t1, t2, t3, t4, t5, gamma)
    factor = p1 - p2 * (1.0 - GQH_GGmax)**p3
    D = damp2[0]/100.0 - factor[0] * GQH_damp_org[0]
    res = least_squares(resDamp,initProps, args=(t1, t2, t3, t4, t5, tau_max, G_max, gamma2, damp2/100.0 - D), bounds=([0.0, 0.0, 0.0],[np.inf, np.inf, np.inf]))
    p1 = res.x[0]
    p2 = res.x[1]
    p3 = res.x[2]

    return t1, t2, t3, t4, t5, p1, p2, p3, D

# calibration function for GQH model based on Vucetic curves - using least squares method
def calibrateGQH_VuceticDobry_LS(PI, tau_max, G_max):
    # get the reference modulus reduction curve
    GGmax2, damp2, gamma2 = VuceticDobry(PI)

    # set a cut-off strain for optimization
    eps_cutoff = 100.0
    index = np.where(gamma2 > eps_cutoff)
    GGmax = np.delete(GGmax2, index)
    gamma = np.delete(gamma2, index)

    # calculate the theta_tau for the reference GGmax curve
    gr = tau_max/G_max
    Gam = (gamma/gr)
    tt = (GGmax + GGmax * Gam - 1.0) / (GGmax**2.0 * Gam)
    tt[np.where(tt > 1.0)] = 1.0

    # initial guess
    t10 = -1.0
    t20 = -3.0
    t30 = 1.0
    t50 = 0.5

    # perform the least squares analysis
    initProp = [t10,t20,t30,t50]
    # res = least_squares(GHQ_tt_res, initProp, args=(tt, tau_max, G_max, gamma), bounds=([-20.0, -20.0, 0.0001, 0.0],[20.0, 20.0, np.inf, 0.99]))
    # t1 = res.x[0]
    # t2 = res.x[1]
    # t3 = res.x[2]
    # t4 = 1.0
    # t5 = res.x[3]
    res = fmin_slsqp(errorGGmax, initProp, args=(tau_max, G_max, gamma2, GGmax2), bounds=[(-20,20),(-20,20),(0,np.inf),(0.1,0.99)], ieqcons=[constraint1, constraint2, constraint3], iprint=-1)
    t1 = res[0]
    t2 = res[1]
    t3 = res[2]
    t4 = 1.0
    t5 = res[3]

    # check the strength and the optimized values
    GQH_GGmax = GQH_modulusreduction(tau_max,G_max, t1, t2, t3, t4, t5, gamma2)
    # app_strength = (G_max*GQH_GGmax*gamma2)[-1]
    # if  (app_strength < 0.95 * tau_max) or (t1 + t2 > 1.0):
    #     initProp = [t10, t20]
    #     res = least_squares(GHQ_tt_res95_1, initProp, args=(tt, tau_max, G_max, gamma, t3, t5), bounds=([-20.0, -20.0],[20.0, 20.0]))
    #     t1 = res.x[0]
    #     t2 = res.x[1]
    #     initProp = [t30, t50]
    #     res = least_squares(GHQ_tt_res95_2, initProp, args=(tt, tau_max, G_max, gamma, t1, t2), bounds=([0.0001, 0.0],[np.inf, 0.99]))
    #     t3 = res.x[0]
    #     t5 = res.x[1]
    
    # optimize the damping parameters
    initProps = [1.0,0.0,1.0]
    res = least_squares(resDamp,initProps, args=(t1, t2, t3, t4, t5, tau_max, G_max, gamma2, damp2/100.0), bounds=([0.0, 0.0, 0.0],[np.inf, np.inf, np.inf]))
    p1 = res.x[0]
    p2 = res.x[1]
    p3 = res.x[2]

    # calculate the minimum damping and revise the parameters
    GQH_damp_org  = GQH_damping(tau_max, G_max, t1, t2, t3, t4, t5, gamma)
    factor = p1 - p2 * (1.0 - GQH_GGmax)**p3
    D = damp2[0]/100.0 - factor[0] * GQH_damp_org[0]
    res = least_squares(resDamp,initProps, args=(t1, t2, t3, t4, t5, tau_max, G_max, gamma2, damp2/100.0 - D), bounds=([0.0, 0.0, 0.0],[np.inf, np.inf, np.inf]))
    p1 = res.x[0]
    p2 = res.x[1]
    p3 = res.x[2]

    return t1, t2, t3, t4, t5, p1, p2, p3, D

# calibration function for GQH model based on Zhang 2005 curves - using least squares method
def calibrateGQH_Zhang2005_LS(p, PI, tau_max, G_max, age='Quaternary', Patm=101.3):
    # get the reference modulus reduction curve
    GGmax2, damp2, gamma2 = Zhang2005(p, PI, age, Patm)

    # set a cut-off strain for optimization
    eps_cutoff = 100.0
    index = np.where(gamma2 > eps_cutoff)
    GGmax = np.delete(GGmax2, index)
    gamma = np.delete(gamma2, index)

    # calculate the theta_tau for the reference GGmax curve
    gr = tau_max/G_max
    Gam = (gamma/gr)
    tt = (GGmax + GGmax * Gam - 1.0) / (GGmax**2.0 * Gam)
    tt[np.where(tt > 1.0)] = 1.0

    # initial guess
    t10 = -1.0
    t20 = -3.0
    t30 = 1.0
    t50 = 0.5

    # perform the least squares analysis
    initProp = [t10,t20,t30,t50]
    # res = least_squares(GHQ_tt_res, initProp, args=(tt, tau_max, G_max, gamma), bounds=([-20.0, -20.0, 0.0001, 0.0],[20.0, 20.0, np.inf, 0.99]))
    # t1 = res.x[0]
    # t2 = res.x[1]
    # t3 = res.x[2]
    # t4 = 1.0
    # t5 = res.x[3]
    res = fmin_slsqp(errorGGmax, initProp, args=(tau_max, G_max, gamma2, GGmax2), bounds=[(-20,20),(-20,20),(0,np.inf),(0.1,0.99)], ieqcons=[constraint1, constraint2, constraint3], iprint=-1)
    t1 = res[0]
    t2 = res[1]
    t3 = res[2]
    t4 = 1.0
    t5 = res[3]

    # check the strength and the optimized values
    GQH_GGmax = GQH_modulusreduction(tau_max,G_max, t1, t2, t3, t4, t5, gamma2)
    # app_strength = (G_max*GQH_GGmax*gamma2)[-1]
    # if  (app_strength < 0.95 * tau_max) or (t1 + t2 > 1.0):
    #     initProp = [t10, t20]
    #     res = least_squares(GHQ_tt_res95_1, initProp, args=(tt, tau_max, G_max, gamma, t3, t5), bounds=([-20.0, -20.0],[20.0, 20.0]))
    #     t1 = res.x[0]
    #     t2 = res.x[1]
    #     initProp = [t30, t50]
    #     res = least_squares(GHQ_tt_res95_2, initProp, args=(tt, tau_max, G_max, gamma, t1, t2), bounds=([0.0001, 0.0],[np.inf, 0.99]))
    #     t3 = res.x[0]
    #     t5 = res.x[1]
    
    # optimize the damping parameters
    initProps = [1.0,0.0,1.0]
    res = least_squares(resDamp,initProps, args=(t1, t2, t3, t4, t5, tau_max, G_max, gamma2, damp2/100.0), bounds=([0.0, 0.0, 0.0],[np.inf, np.inf, np.inf]))
    p1 = res.x[0]
    p2 = res.x[1]
    p3 = res.x[2]

    # calculate the minimum damping and revise the parameters
    GQH_damp_org  = GQH_damping(tau_max, G_max, t1, t2, t3, t4, t5, gamma)
    factor = p1 - p2 * (1.0 - GQH_GGmax)**p3
    D = damp2[0]/100.0 - factor[0] * GQH_damp_org[0]
    res = least_squares(resDamp,initProps, args=(t1, t2, t3, t4, t5, tau_max, G_max, gamma2, damp2/100.0 - D), bounds=([0.0, 0.0, 0.0],[np.inf, np.inf, np.inf]))
    p1 = res.x[0]
    p2 = res.x[1]
    p3 = res.x[2]

    return t1, t2, t3, t4, t5, p1, p2, p3, D
    

def main():
    PI = 50.0
    G_max = 50000.0
    tau_max = 100.0
    # t1,t2,t3,t4,t5,p1,p2,p3,D = calibrateGQH_Darendeli(p,PI,OCR,tau_max,G_max, Patm=2117.0)
    # GGmax, damp, gamma = VuceticDobry(PI)
    # t1,t2,t3,t4,t5,p1,p2,p3,D = calibrateGQH_VuceticDobry_LS(PI,tau_max,G_max)
    GGmax, damp, gamma = Zhang2005(200.0, PI, Patm=101.3)
    t1,t2,t3,t4,t5,p1,p2,p3,D = calibrateGQH_Zhang2005_LS(200.0, PI,tau_max,G_max, Patm=101.3)

    GQH_GGmax = GQH_modulusreduction(tau_max,G_max, t1, t2, t3, t4, t5, gamma)
    GQH_damp_org  = GQH_damping(tau_max, G_max, t1, t2, t3, t4, t5, gamma)
    factor = p1 - p2 * (1.0 - GQH_GGmax)**p3
    GQH_damp = D + factor * GQH_damp_org

    print(" t1 = {0} \n t2 = {1} \n t3 = {2} \n t4 = {3} \n t5 = {4} \n p1 = {5} \n p2 = {6} \n p3 = {7} \n D = {8}".format( \
        t1, t2, t3, t4, t5, p1, p2, p3, damp[0]-GQH_damp[0]))

    plt.figure("GGmax")
    plt.semilogx(gamma*100.0, GGmax, label='Vucetic')
    plt.semilogx(gamma*100.0, GQH_GGmax, label='GQH')
    plt.xlabel('gamma (%)')
    plt.ylabel('G/Gmax')
    plt.legend()
    plt.grid(which='both')
    plt.show(block=False)
    
    plt.figure("Stress-Strain")
    plt.plot(gamma*100.0, GGmax*G_max*gamma, label='Vucetic')
    plt.plot(gamma*100.0, GQH_GGmax*G_max*gamma, label='GQH')
    plt.plot(gamma*100.0, 0.0*gamma+tau_max, 'k--')
    plt.plot(gamma*100.0, 0.0*gamma+(0.95*tau_max), 'k--')
    plt.xlabel('gamma (%)')
    plt.ylabel('tau (psf)')
    plt.legend()
    plt.grid(which='both')
    plt.show(block = False)

    plt.figure("Damp")
    plt.semilogx(gamma*100.0, damp, label='Vucetic')
    plt.semilogx(gamma*100.0, GQH_damp_org*100.0, label='GQH_original')
    plt.semilogx(gamma*100.0, GQH_damp*100.0, label='GQH_corrected')
    plt.ylabel('damping ratio (%)')
    plt.xlabel('gamma (%)')
    plt.legend()
    plt.grid(which='both')
    plt.show()

    p = 650.0
    PI = 0.0
    OCR = 1.0
    G_max = 130/32.2 * 1785.21**2.0
    tau_max = 20000.0
    GGmax, damp, gamma = Darendeli(p, PI, OCR, Patm=2117.0)
    # t1,t2,t3,t4,t5,p1,p2,p3,D = calibrateGQH_Darendeli(p,PI,OCR,tau_max,G_max, Patm=2117.0)
    t1,t2,t3,t4,t5,p1,p2,p3,D = calibrateGQH_Darendeli_LS(p,PI,OCR,tau_max,G_max, Patm=2117.0)

    GQH_GGmax = GQH_modulusreduction(tau_max,G_max, t1, t2, t3, t4, t5, gamma)
    GQH_damp_org  = GQH_damping(tau_max, G_max, t1, t2, t3, t4, t5, gamma)
    factor = p1 - p2 * (1.0 - GQH_GGmax)**p3
    GQH_damp = D + factor * GQH_damp_org

    print(" t1 = {0} \n t2 = {1} \n t3 = {2} \n t4 = {3} \n t5 = {4} \n p1 = {5} \n p2 = {6} \n p3 = {7} \n D = {8}".format( \
        t1, t2, t3, t4, t5, p1, p2, p3, damp[0]-GQH_damp[0]))

    plt.figure("GGmax")
    plt.semilogx(gamma*100.0, GGmax, label='Darendeli')
    plt.semilogx(gamma*100.0, GQH_GGmax, label='GQH')
    plt.xlabel('gamma (%)')
    plt.ylabel('G/Gmax')
    plt.legend()
    plt.grid(which='both')
    plt.show(block=False)
    
    plt.figure("Stress-Strain")
    plt.plot(gamma*100.0, GGmax*G_max*gamma, label='Darendeli')
    plt.plot(gamma*100.0, GQH_GGmax*G_max*gamma, label='GQH')
    plt.plot(gamma*100.0, 0.0*gamma+tau_max, 'k--')
    plt.plot(gamma*100.0, 0.0*gamma+(0.95*tau_max), 'k--')
    plt.xlabel('gamma (%)')
    plt.ylabel('tau (psf)')
    plt.legend()
    plt.grid(which='both')
    plt.show(block = False)

    plt.figure("Damp")
    plt.semilogx(gamma*100.0, damp, label='Darendeli')
    plt.semilogx(gamma*100.0, GQH_damp_org*100.0, label='GQH_original')
    plt.semilogx(gamma*100.0, GQH_damp*100.0, label='GQH_corrected')
    plt.ylabel('damping ratio (%)')
    plt.xlabel('gamma (%)')
    plt.legend()
    plt.grid(which='both')
    plt.show()

    p = 650.0
    Cu = 2.0
    D50 = 0.1
    G_max = 130/32.2 * 1000.0**2.0
    tau_max = 10000.0
    GGmax, damp, gamma = Menq(p, Cu, D50, Patm=2117.0)
    t1,t2,t3,t4,t5,p1,p2,p3,D = calibrateGQH_Menq_LS(p, Cu, D50,tau_max,G_max, Patm=2117.0)

    GQH_GGmax = GQH_modulusreduction(tau_max,G_max, t1, t2, t3, t4, t5, gamma)
    GQH_damp_org  = GQH_damping(tau_max, G_max, t1, t2, t3, t4, t5, gamma)
    factor = p1 - p2 * (1.0 - GQH_GGmax)**p3
    GQH_damp = D + factor * GQH_damp_org

    print(" t1 = {0} \n t2 = {1} \n t3 = {2} \n t4 = {3} \n t5 = {4} \n p1 = {5} \n p2 = {6} \n p3 = {7} \n D = {8}".format( \
        t1, t2, t3, t4, t5, p1, p2, p3, damp[0]-GQH_damp[0]))

    plt.figure("GGmax")
    plt.semilogx(gamma*100.0, GGmax, label='Menq')
    plt.semilogx(gamma*100.0, GQH_GGmax, label='GQH')
    plt.xlabel('gamma (%)')
    plt.ylabel('G/Gmax')
    plt.legend()
    plt.grid(which='both')
    plt.show(block=False)
    
    plt.figure("Stress")
    plt.plot(gamma*100.0, GGmax*G_max*gamma, label='Menq')
    plt.plot(gamma*100.0, GQH_GGmax*G_max*gamma, label='GQH')
    plt.plot(gamma*100.0, 0.0*gamma+tau_max, 'k--')
    plt.plot(gamma*100.0, 0.0*gamma+(0.95*tau_max), 'k--')
    plt.xlabel('gamma (%)')
    plt.ylabel('tau (psf)')
    plt.legend()
    plt.grid(which='both')
    plt.show(block = False)

    plt.figure("Damp")
    plt.semilogx(gamma*100.0, damp, label='Menq')
    plt.semilogx(gamma*100.0, GQH_damp_org*100.0, label='GQH_original')
    plt.semilogx(gamma*100.0, GQH_damp*100.0, label='GQH_corrected')
    plt.ylabel('damping ratio (%)')
    plt.xlabel('gamma (%)')
    plt.legend()
    plt.grid(which='both')
    plt.show()


if __name__ == '__main__':
    main()
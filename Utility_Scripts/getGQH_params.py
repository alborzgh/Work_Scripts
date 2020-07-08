# get GQH parameters
from calibrateGQH import calibrateGQH_Darendeli

PI = 0.0
OCR = 1.0
p = 650.0
G_max = 130/32.2*1000**2.0
tau_max = 1000.0

t1,t2,t3,t4,t5,p1,p2,p3,D = calibrateGQH_Darendeli(p,PI,OCR,tau_max,G_max, Patm=2117.0)
print("t1 = {0} \n t2 = {1} \n t3 = {2} \n t4 = {3} \n t5 = {4} \n p1 = {5} \n p2 = {6} \n p3 = {7}".format( \
        t1, t2, t3, t4, t5, p1, p2, p3))
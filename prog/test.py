import numpy as np
import matplotlib.pyplot as plt

q = 1.6*pow(10,-19)
B_0 = 40/10000
m = 9.1*pow(10,-31)
r_0 = 1.0
W = 1.6*pow(10,-14)
r_s = 1.065

nu1 = pow(r_0,2)*2*W/m
nu2 = q*B_0*r_0/m

C1 = nu1/pow(r_s,3)
C2 = nu2*nu2*np.log(r_s/r_0)
C3 = (nu1*3/r_s - nu2*nu2*(1-np.log(r_s/r_0)))/pow(r_s,2)
C = C1 - C2 + C3*r_s
K2 = r_0 - C3/C

r=(r_0 - C3/C)*np.cos(np.sqrt(C3)*t) + C3/C

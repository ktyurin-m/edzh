import numpy as np
import matplotlib.pyplot as plt

q = 1.6*pow(10,-19)
B_0 = 40/10000
m = 9.1*pow(10,-31)
r_0 = 1.0
W = 1.6*pow(10,-14)

nu1 = pow(r_0,2)*2*W/m
nu2 = q*B_0*r_0/m
z_v = 0

def funcz(r):
    return  0+ nu2*np.log(r/r_0)
def funcr(r):
    return (nu1/(pow(r,3.0)) - nu2*funcz(r)/r)
r_s = 1.0662 

K = np.sqrt(3*nu1/pow(r_s,2) + pow(nu2,2)*(1-np.log(r_s/r_0)))/r_s

C3 = nu1/pow(r_s,3) - pow(nu2,2)*np.log(r_s/r_0)/r_s
C2 = r_0 - C3/pow(K,2) - r_s

def an(t):
    return r_s + C2*np.cos(K*t) + C3/pow(K,2)
def tr(t):
    return np.sin(K*t)*C2/K + C3*t/pow(K,2)
def fz(t):
    return nu2*np.log(r_s/r_0)*t + tr(t)*nu2/r_s

def theta(r):
    return np.sqrt(2*W/m)*r_0/pow(r,2)

def vcd(r):
    return 2*W*r_0/(q*B_0*pow(r,2))

def anvz(t):
    return nu2*(np.log(r_s/r_0)+(C2*np.cos(K*t) + C3/pow(K,2))/r_s)

R_t = []
t = []
ar3 = []
anf = []
ar5 = []
ar6x = []
ar6y = []
ar7x = []
ar7y = []
ar8v = []
ar8v1 = []

arz1 = []
arz2 = []
flg,axs = plt.subplots(nrows= 1, ncols= 4)

print(
    "z(velocity)= "+ str(2*W/(q*B_0))
    )

def ur():
    v = 0
    r_0 = 1
    T = m*2*np.pi/(q*B_0)
    n = 10
    dt =  T/10000
    z = 0
    v_z = 0
    thet = 0
    R_t.append(r_0)
    t.append(0)
    ar3.append(z)
    anf.append(an(0))
    ar5.append(fz(0))
    ar6x.append(r_0)
    ar6y.append(0)
    
    ar7x.append(r_0)
    ar7y.append(0)
    ar8v.append(vcd(r_0))
    ar8v1.append(0)
    arz1.append(0)
    arz2.append(0)
    
    vzz = 0
    vv1 = 0
    i = 0
    while(True):
        v_z = funcz(r_0)
        z = z + funcz(r_0)*dt

        v = v + dt*funcr(r_0)
        r_0 = r_0 + v*dt
        
        thet = thet + theta(r_0)*dt
        vv1 = vv1 + vcd(r_0)*dt
        ar6x.append(r_0*np.cos(thet))
        ar6y.append(r_0*np.sin(thet))
        if(thet >= np.pi * 2):
            axs[0].plot(t,R_t,'g',t,anf,'r')
            #,t,anf,'r'
            axs[0].grid()
            axs[0].set_xlabel("$t$, с")
            axs[0].set_ylabel("$r(t)$, м", loc="top")

            axs[1].plot(t,ar3,'g',t,ar5,'r',t,ar8v1,'b')
            #,t,ar5,'r', t, ar8v, 'b', t, ar8v1,'y'
            axs[1].grid()
            axs[1].set_xlabel("$t$, c")
            axs[1].set_ylabel("$z(t)$, м", loc="top")
            
            axs[2].plot(ar6x,ar6y)
            axs[2].set_xlabel("$x$")
            axs[2].set_ylabel("$y$")
            axs[2].grid()
            # axs[2].plot(t,anf)

            axs[3].plot(t, arz1,'g', t, arz2,'r')
            axs[3].set_xlabel("$t$, c")
            axs[3].set_ylabel("$v_z(t)$, м", loc="top")
            axs[3].grid()
            print(i)
            plt.show()
            R_t.clear()
            t.clear()
            break
        i = i + dt
        ar3.append(z)
        ar8v.append(vcd(r_0))
        ar8v1.append(vv1)
        R_t.append(r_0)
        t.append(i)
        anf.append(an(i))
        ar5.append(fz(i))
        arz1.append(v_z)
        arz2.append(anvz(i))
ur()
#code can be found at:
# https://github.com/mondalbhaskar/AE-502-Homework-Project-3.git

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import ode
from scipy.integrate import LSODA
from scipy.integrate import solve_ivp

def eccentricanomalyFROMmeananomaly(M,e):
    #Get eccentric anomaly using iteration
    #set mean anomaly as first guess
    E = M
    tolerance =1e-8
    for i in range(1000):
        f = E -e*math.sin(E)-M
        fprime = 1-e*math.cos(E)
        dE = -f/fprime
        E = E+dE
        
        if abs(dE) < tolerance:
            #print("Obtained Ecentric anomaly from Mean anomaly")
            break
    #now we have eccentric anomaly
    return E

def convertToCartesian(a,e,i,w,W,M):
    E = eccentricanomalyFROMmeananomaly(M, e)
    theta = 2*math.atan2( (1+e)**0.5 * math.sin(E/2),\
                                (1-e)**0.5 * math.cos(E/2))
    h = math.sqrt(mu*a*(1-e**2))
    r = (h**2/mu)/(1-e*math.cos(theta))

    x = r*(math.cos(W)*math.cos(w+theta)-math.sin(W)*math.sin(w+theta)*math.cos(i))
    y = r*(math.sin(W)*math.cos(w+theta)+math.cos(W)*math.sin(w+theta)*math.cos(i))
    z = r*(math.sin(i)*math.sin(w+theta))
    return x,y,z
    

omega_rot = 0.01 # frame rotation
a = 1 # semi-major axis
e = 0.5 # eccentricity
i = 45*math.pi/180 # inclination
mu = 1
#Time derivatives of the orbital elements in Delauney action angle variables
# L,G,H are the actions
#l,g,h are the angles
def dL_dt(L,G,H,l,g,h):
    return 0
def dG_dt(L,G,H,l,g,h):
    return 0
def dH_dt(L,G,H,l,g,h):
    return 0
def dl_dt(L,G,H,l,g,h):
    return 1/(L**3)
def dg_dt(L,G,H,l,g,h):
    #return 1/(2*G**2)
    return omega_rot*math.cos(i)
def dh_dt(L,G,H,l,g,h):
    return omega_rot

def f(t, y):
    dydt = np.zeros(6)
    # y[0] = L, y[1] = G, y[2] = H, y[3] = l, y[4] = g, y[5] = h
    dydt[0] = dL_dt(y[0],y[1],y[2],y[3],y[4],y[5])
    dydt[1] = dG_dt(y[0],y[1],y[2],y[3],y[4],y[5])
    dydt[2] = dH_dt(y[0],y[1],y[2],y[3],y[4],y[5])
    dydt[3] = dl_dt(y[0],y[1],y[2],y[3],y[4],y[5])
    dydt[4] = dg_dt(y[0],y[1],y[2],y[3],y[4],y[5])
    dydt[5] = dg_dt(y[0],y[1],y[2],y[3],y[4],y[5])
    return dydt

totalTime=100 #time units



# Set the initial conditions
L=[1.0]
G=[L[0]*math.sqrt(1-e**2)]
H=[G[0]*math.cos(i)]
l=[0] #mean anomaly
g=[0] #longitude of ascending node
h=[0] #argument of pericenter
#time=[0] #array for time

t=0
#solve ivp
# # Create the `ode` object with `lsoda` solver
# integrator = ode(f).set_integrator('lsoda', rtol=1e-6, atol=1e-9, max_step=0.0010)

# integrator.set_initial_value([L[0],G[0],H[0],l[0],g[0],h[0]], 0)

# #integrate the system
# while integrator.successful() and integrator.t < totalTime:
#     integrator.integrate(integrator.t+1)
#     L.append(integrator.y[0])
#     G.append(integrator.y[1])
#     H.append(integrator.y[2])
#     l.append(integrator.y[3])
#     g.append(integrator.y[4])
#     h.append(integrator.y[5])
#     time.append(integrator.t)
#     t+=1

#dt = 0.0001
time = np.linspace(0,totalTime,10000,endpoint=True)
dt=time[2]-time[1]
#Propagate the system in time
for ttime in time:
    L.append(L[t]+dL_dt(L[t],G[t],H[t],l[t],g[t],h[t])*dt)
    G.append(G[t]+dG_dt(L[t],G[t],H[t],l[t],g[t],h[t])*dt)
    H.append(H[t]+dH_dt(L[t],G[t],H[t],l[t],g[t],h[t])*dt)
    l.append(l[t]+dl_dt(L[t],G[t],H[t],l[t],g[t],h[t])*dt)
    g.append(g[t]+dg_dt(L[t],G[t],H[t],l[t],g[t],h[t])*dt)
    h.append(h[t]+dh_dt(L[t],G[t],H[t],l[t],g[t],h[t])*dt)
    #time.append(time[t]+dt)
    t+=1



#Convert to Keplerian system
n=math.sqrt(mu/np.power(a,3))
a_kep = np.empty(len(L),dtype=float)
e_kep = np.empty(len(L),dtype=float)
i_kep = np.empty(len(L),dtype=float)
Omega_kep = np.empty(len(L),dtype=float)

M = l # mean anomaly
omega_kep = g
for index in range(len(L)):
    a_kep[index] = np.square(L[index])
    e_kep[index] = np.sqrt(1-(np.square(G[index]))/(np.square(L[index])))
    i_kep[index] = np.arccos(H[index]/(G[index]))
    Omega_kep[index] = h[index]-g[index]

X = np.empty(len(L),dtype=float)
Y = np.empty(len(L),dtype=float)
Z = np.empty(len(L),dtype=float)

for index in range(len(L)):
    X[index],Y[index],Z[index] = convertToCartesian(a_kep[index],e_kep[index],i_kep[index],omega_kep[index],Omega_kep[index],M[index])




#Plot the orbit
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(X, Y, Z)

#show point (0,0,0) as a red dot, which is the central body
ax.scatter(0,0,0,c='r',marker='o')

#show the initial position as a green dot
ax.scatter(X[0],Y[0],Z[0],c='g',marker='o')
#show the final position as a blue dot
ax.scatter(X[-1],Y[-1],Z[-1],c='b',marker='o')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')


#save the figure as a jpg
fig.savefig('orbit2.jpg')

plt.show()





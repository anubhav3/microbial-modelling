"""
Title: Modelling microbial communities
Author: Anubhav Gupta
Date: 31/10/2018
"""


from random import random,seed
from scipy.integrate import ode
import numpy as np
import matplotlib.pyplot as plt


n = 5 #Number of species
n_population = 10000. #Total number of initial populations
n_resource = 300. #Number of initial resource
mu_max = [0.07 for i in range(n)] #Maximum specific growth rate
half_sat = [5. for i in range(n)] #Half Saturation Constant
yieldd = 10**7. #Number of cells per unit volume of a resource
m = [0.03 for i in range(n)] #Mortality Rate


#seed(3) #Setting the seed for generating seqeunce of random numbers


#Model for the population
def model(N):
	dNdt = [0 for i in range(n+1)]
	for i in range(n):
		dNdt[i] = mu_max[i]*N[n]/(N[n]+half_sat[i])*N[i] - m[i]*N[i]
	dNdt[n] = -(1/yieldd)*mu_max[1]*N[n]/(N[n]+half_sat[1])*sum(N[:n])	#This line would not be the same if the mu_max and half_sar are not same for all species.
	return dNdt


#The Runge-Kutta 4th Order Integration method (RK4)
def rk(h,N):
	k1 = np.multiply(model(N),h)
	k2 = np.multiply(model(N+k1/2),h)
	k3 = np.multiply(model(N+k2/2),h)
	k4 = np.multiply(model(N+k3),h)
	yn = N + (1./6)*(k1+2*k2+2*k3+k4)
	if(yn[n]<=0):  #Just to ensure no. of resources is always non-negative
		yn[n] = 0
	return yn

#Initial Populations
IN = [10**(2+5*random()) for i in range(n)]
#Initial resource (Using a single array to store populations and resource)
IN.append(n_resource)


h = 0.05 #Integration step size
t_final = 400
ts =[]
ts.append(IN)
N = IN
t = 0.05
while(t<t_final):
	N = rk(h,N)   #Applying RK4 method
	ts.append(N)
	t+=h
ts = np.transpose(ts)

x_t = [t for t in range(int(len(ts[1])))]

#Just to check the final population
for i in range(n+1):
	print(ts[i][len(x_t)-1])



#Plotting
fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('time (t)')
ax1.set_ylabel('populations', color=color)
for i in range(n):	
	plt.plot(x_t,ts[i])

ax2 = ax1.twinx()
color = 'tab:blue'
ax2.set_ylabel('resource', color="black")	
ax2.plot(x_t,ts[n],color="black")
# fig.tight_layout()
plt.show()






#!/usr/bin/python
# -*- coding: utf8 -*-

import numpy as np


def rk4(f,dt,t,y):
  k1 = f(t, y)
  k2 = f(t + dt/2., y + (dt/2.)*k1)
  k3 = f(t + dt/2., y + (dt/2.)*k2)
  k4 = f(t + dt   , y +  dt    *k3)
  return (dt*(k1+2.*k2+2.*k3+k4)/6.)

def f(t, y):
  def f1(t, y):
    return (y[2])
  def f2(t, y):
#    return (mu * (1.-y[1]**2) * y[2] - w**2 * y[1] + a * np.sin(Wf * t))
    return (-c*y[2] - np.sin(y[1]) + a*np.cos(Wf*t)) 
  return np.array([0.,f1(t,y),f2(t,y)])
# osc arm  (-(w**2)*np.sin(y1)-c*y2+a*np.sin(Wf*t))
# van der pol (mu*(1.-y1**2)*y2-w**2*y1+a*np.sin(Wf*t))

#################
# condicion inicial y definicion del vector (lista) con las variables

# uso un placeholder para que y1 sea y[1], y2=y[2] etc La coordenada 0 de y no tiene valor
y=np.array([0.,1.,0.])
# el tiempo es explicito
t=0;

dt=np.pi/300.
tf=10000.
m=int((tf-t)/dt+0.5)+1
nc=100000.
## parámetros de la edo
# w=1.
# mu=5.
# a=9.1
# Wf=(9.53)
# fv=open("vanderpol_{0:05d}_{1:04d}_{2:04d}_{3:02d}{4:02d}.dat".format(int(mu*1000),int(Wf*1000),int(a*1000),int(y[1]*10),int(y[2]*10)),"w")

# ## oscilador-forzado

w = 1.
c = 0.25
Wf= 2./3.
delta=dt/100000.
goal=3.*np.pi
a=0.
fv=open("pendulo_{0:02d}_{1:04d}_{2:04d}.bif".format(int(c*100),int(Wf*1000),int(0*1000)),"w")
fv.close()

pi2=2.*np.pi


for a in np.arange(0.5,1.00,0.001):
  fv=open("pendulo_{0:02d}_{1:04d}_{2:04d}.bif".format(int(c*100),int(Wf*1000),int(1*1000)),"a")
  t=0
  i=0
  n=1
  y=np.array([0.,1.,0.])
  while (n<151):
    i += 1
    if (abs(t-n*goal)<delta):
      print a,n
      n+=1
      if (n>100):
          fv.write("{0:.03f} {1:.2f} {2:.7f} {3:.7f}\n".format(a,t,y[1], y[2]))
    y = y+rk4(f,dt,t,y)
    t += dt
    if y[1] > np.pi:
      y[1] -= pi2
    else:
      if y[1] < -np.pi:
        y[1] += pi2

  fv.close()

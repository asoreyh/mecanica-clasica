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
    return (-c*y[2] - w**2 * np.sin(y[1]) + a*np.cos(Wf*t)) 
  return np.array([0.,f1(t,y),f2(t,y)])
# osc arm  (-(w**2)*np.sin(y1)-c*y2+a*np.sin(Wf*t))
# van der pol (mu*(1.-y1**2)*y2-w**2*y1+a*np.sin(Wf*t))

#################
# condicion inicial y definicion del vector (lista) con las variables

# uso un placeholder para que y1 sea y[1], y2=y[2] etc La coordenada 0 de y no tiene valor
y=np.array([0.,1.,0.])
# el tiempo es explicito
t=0;

dt=0.005
tf=200.
m=int((tf-t)/dt+0.5)+1

# ## parámetros de la edo
# w=1.
# mu=5.
# a=9.1
# Wf=(9.53)
# fv=open("vanderpol_{0:05d}_{1:04d}_{2:04d}_{3:02d}{4:02d}.dat".format(int(mu*1000),int(Wf*1000),int(a*1000),int(y[1]*10),int(y[2]*10)),"w")

## oscilador-forzado
w = 1.
c = 0.1
Wf= 1.0 # np.pi
a = 0.3
fv=open("pendulo_{0:02d}_{1:04d}_{2:04d}.dat".format(int(c*100),int(Wf*1000),int(a*1000)),"w")
for i in range(m):
  fv.write("{0:d} {1:.2f} {2:.5f} {3:.5f}\n".format(i,t,y[1], y[2]))
  y = y+rk4(f,dt,t,y)
  t += dt
fv.close()

#!/usr/bin/python
# -*- coding: utf8 -*-

import numpy as np

def rk4(f,dt,t,y):
  k1 = f(t, y)
  k2 = f(t + dt/2., y + (dt/2.)*k1)
  k3 = f(t + dt/2., y + (dt/2.)*k2)
  k4 = f(t + dt   , y +  dt    *k3)
  return (dt*(k1+2.*k2+2.*k3+k4)/6.)

## HH no depende explícitamente del tiempo, pero seguimos
def f(t, y):
  def f1(t, y):
    return (y[2])
  def f2(t, y):
    return (-y[1]-l*(y[3]**2-y[1]**2))
  def f3(t, y):
    return (y[4])
  def f4(t, y):
    return (-y[3]-2.*l*y[3]*y[1])
  return np.array([0.,f1(t,y),f2(t,y),f3(t,y),f4(t,y)])

#################
# condicion inicial y definicion del vector (lista) con las variables

# uso un placeholder para que y1 sea y[1], y2=y[2] etc La coordenada 0 de y no tiene valor
# el tiempo es explicito
t=0;
dt=0.005
tf=400.
m=int((tf-t)/dt+0.5)+1

## henon-heiles
l=1.
y=np.array([0.0, 0.020, -0.080, 0.000, 0.000])
y[4]=np.sqrt(2.*1/12.-y[2]**2-y[1]**2+2/3.*y[1]**3)

fv=open("hh_{0:04d}_{1:04d}_{2:04d}.dat".format(int(l*1000),int(y[1]*1000), int(y[3]*1000.)),"w")

for i in range(m):
  fv.write("{0:d} {1:.2f} {2:.9f} {3:.9f} {4:.9f} {5:.9f}\n".format(i, t, y[1], y[2], y[3], y[4]))
  y = y+rk4(f,dt,t,y)
  t += dt
fv.close()

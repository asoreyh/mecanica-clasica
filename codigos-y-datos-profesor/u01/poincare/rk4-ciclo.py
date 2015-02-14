#!/usr/bin/python
# -*- coding: utf8 -*-

import numpy as np
import sys

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
    return (-c*y[2] - w**2 * np.sin(y[1]) + a*np.sin(Wf*t)) 
  return np.array([0.,f1(t,y),f2(t,y)])
# osc arm  (-(w**2)*np.sin(y1)-c*y2+a*np.sin(Wf*t))
# van der pol (mu*(1.-y1**2)*y2-w**2*y1+a*np.sin(Wf*t))

#################
# condicion inicial y definicion del vector (lista) con las variables
# uso un placeholder para que y1 sea y[1], y2=y[2] etc La coordenada 0 de y no tiene valor
y=np.array([0.,0.5,0.])
# el tiempo es explicito
t=0;

# ## parámetros de la edo
# w=1.
# mu=5.
# a=9.1
# Wf=(9.53)
# fv=open("vanderpol_{0:05d}_{1:04d}_{2:04d}_{3:02d}{4:02d}.dat".format(int(mu*1000),int(Wf*1000),int(a*1000),int(y[1]*10),int(y[2]*10)),"w")

## oscilador-forzado
w = 1.
c = 0.25
Wf= 2./3.
a = 1.00
fv=open("pendulo_{0:02d}_{1:04d}_{2:04d}.dat".format(int(c*100),int(Wf*1000),int(a*1000)),"w")
# ciclo del forzado
res=300.
nmax=50


dt=np.pi/res
i=0
n=0
iciclo=2.*res/Wf
# iciclo debe ser entero ( y no aproximadamente entero!)
if not iciclo.is_integer():
  print "Modificar los valores de dt para tener índice de ciclos (iciclo) entero"
  sys.exit()

pi2=2.*np.pi

ipoinc=True

if ipoinc:
  fp=open("pendulo_{0:02d}_{1:04d}_{2:04d}.poi".format(int(c*100),int(Wf*1000),int(a*1000)),"w")
  print "pendulo_{0:02d}_{1:04d}_{2:04d}.poi".format(int(c*100),int(Wf*1000),int(a*1000)),nmax
# ciclo
while (n<nmax):
  if n>10:
    fv.write("{0:d} {1:.2f} {2:.5f} {3:.5f}\n".format(i,t,y[1], y[2]))
  y = y+rk4(f,dt,t,y)
  if y[1] > np.pi:
    y[1] -= pi2
  else:
    if y[1] < -np.pi:
      y[1] += pi2

  i +=1 
  t += dt
  if (i==iciclo):
    print n 
    n += 1
    i=0
    if ipoinc:
      if n>10:
        fp.write("{0:d} {1:.2f} {2:.5f} {3:.5f}\n".format(i,t,y[1], y[2]))
      
fv.close()
if ipoinc:
  fp.close()

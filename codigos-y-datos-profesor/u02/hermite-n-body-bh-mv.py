#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Leapfrog algorithm for the n-body problem
Created on Sun Oct 12 16:21:57 2014
@author: asoreyh
"""

# Numpy!
import numpy as np
from numpy import linalg as la

# Fancy functions

## Gravity, accelerations and ks
def get_ak():
  for i in range(n):
    a[i]=(0.,0.,0.)
    k[i]=(0.,0.,0.)
  mint2 = 1e100
  tji2 = 0.
  for i in range(n):
    for j in range(i+1,n):
      rji=r[j]-r[i]
      rjin=la.norm(rji)
      vji=v[j]-v[i]
      vjin=la.norm(vji)
      aji = G * rji / (rjin**3)
      ajin = la.norm(aji)
      kji = G * (vji / (rjin**3) - 3.* np.dot(rji,vji)* rji / (rjin**5))
      a[i] += m[j]*aji
      a[j] -= m[i]*aji
      k[i] += m[j]*kji
      k[j] -= m[i]*kji
      # time estimation, first step using speeds
      if vjin:
        tji2 = (rjin / vjin)**2
        if (tji2 < mint2):
          mint2 = tji2
      # time estimation, second step using free fall
      tji2 = rjin/ajin
      if (tji2 < mint2):
        mint2 = tji2
  return np.sqrt(mint2)

## Predict
def predict():
  for i in range(n):
    r[i] += v[i] * dt + a[i] * (dt**2)/2. + k[i] * (dt**3)/6.
    v[i] += a[i] * dt + k[i] * (dt**2)/2.

## Correct
def correct():
  for i in range(n):
    v[i] = vp[i] + (ap[i] + a[i]) * dt / 2. + (kp[i] - k[i])*(dt**2)/12.
    r[i] = rp[i] + (vp[i] + v[i]) * dt / 2. + (ap[i] - a[i])*(dt**2)/12.


## Kinetic energy
def e_k():
  e = 0.
  for i in range(n):
    e += 0.5*m[i]*la.norm(v[i])**2
  return e
  
## Potential energy
def e_p():
  e=0.
  for i in range(n):
    for j in range(i+1,n):
      rji=r[j]-r[i]
      e -= G * m[i] * m[j] / la.norm(rji)
  return e

## Total energy
def e():
  return (e_k()+e_p())

## Output
def output(t,r,v):
  for b in range(n):
    print t,
    for c in range(dim):
      print r[b][c],
    print
  print
  print

## Initial Conditions
def initialize():
  r2=R**2
  cs=0
  while (cs<n):
    x = np.random.uniform(-1., 1.)
    y = np.random.uniform(-1., 1.)
    z = np.random.uniform(-1., 1.)
    if (x*x+y*y+z*z<=r2):
      m[cs] = mbase
      if not fixmass:
        m[cs] = (int(cs)%10)+1.
      r[cs] = (x, y, z)
      v[cs] = (0, 0, 0)
      cs += 1

# times
t=0.
tf=1.
dt=0.
ts=0.
tg=0.02
tout=0.
tsteps=20 # we have an output on intervals tg/tsteps

# Constants
#G=6.67384e-11
G=1.
dim=3
R=1.
n=51
step=0

# Holding vectors
m=np.zeros(n)          # mass
# current
r=np.zeros((n,dim))    # position
v=np.zeros((n,dim))    # velocity
a=np.zeros((n,dim))    # acceleration
k=np.zeros((n,dim))    # da/dt=d3r/dt3
# old
rp=np.zeros((n,dim))    # position
vp=np.zeros((n,dim))    # velocity
ap=np.zeros((n,dim))    # acceleration
kp=np.zeros((n,dim))    # da/dt=d3r/dt3

# Initial conditions
# Starting from an spherical configuration
# including here cold-start.py
fixmass=False
mbase=1.
mdev=0.2
initialize()

# now, first start become a blackhole at center
m[0] = 40
r[0] = (0, 0, 0)
v[0] = (0, 0, 0)

# saving initial energy
e_i=e()

# Algorithm
# Here we go. We need to introduce some changes from 
# leapfrog algorithm, mainly because we have to calculate
# not only a but also k=d3r/dt3

# start by calculating a and k
ts = get_ak()
# time loop
while (t<=tf):
  if(t>=tout):
    output(t,r,v)
    tout += tg/tsteps

  for i in range(n):
    rp[i] = r[i]
    vp[i] = v[i]
    ap[i] = a[i]
    kp[i] = k[i]
  dt = tg * ts;
  predict()
  ts = get_ak()
  correct() 
  t += dt
  step += 1

# Good bye
output(t,r,v)
e_f=e()
print "# Energía inicial:", e_i
print "# Energía final:", e_f
print "# Var relativa:", (e_f-e_i)/e_i

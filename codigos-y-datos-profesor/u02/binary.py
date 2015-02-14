#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Predictor-Corrector method based on Hermite Schema
with variable time interval
Calculate n-body interactions up to 4h order
Created on Sun Oct 12 16:21:57 2014
@author: asoreyh
"""

# Numpy
import numpy as np
from numpy import linalg as la

# Fancy functions

## Gravity, accelerations and jerks
def get_ak():
  for i in range(n):
    a[i]=(0.,0.,0.)
    k[i]=(0.,0.,0.)
  mint2 = 1e100
  tji2 = 0.
  for i in range(n):
    for j in range(i+1,n):
      rji=r[j]-r[i]      # relative position ij pair
      rjin=la.norm(rji)  # distance ij
      vji=v[j]-v[i]      # relative velocity 
      vjin=la.norm(vji)  # relative speed
      aji = G * rji / (rjin**3)  # internal acceleration vector
      ajin = la.norm(aji)        # and modulus
      # jerk, d3r/dt3
      kji = G * (vji / (rjin**3) - 3.* np.dot(rji,vji)* rji / (rjin**5))
  
      a[i] += m[j]*aji    # acceleration over body i
      a[j] -= m[i]*aji    # antisymmetric matrix -> over body j
      k[i] += m[j]*kji    # same for jerks
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
  return np.sqrt(mint2) # return with the minimum time

## Predict
## Predictor step: Euler Algorithm using taylor series up to 3th order
def predict(): 
  for i in range(n):
    r[i] += v[i] * dt + a[i] * (dt**2)/2. + k[i] * (dt**3)/6.
    v[i] += a[i] * dt + k[i] * (dt**2)/2.

## Correct
## Corrector step: Hermite Schema up to 4th order
## first calculate velocity, since it is needed for position
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

## orbital temperature

def to():
  l1=1.
  l2=1.
  torb = pow((1./(16*np.pi*sigma))*( l1/(la.norm(r[0]-r[2]))**2 + l2/(la.norm(r[1]-r[2]))**2 ),0.25)
  return torb

## Output
## usefull for gnuplot index method on 'a' steps
## a=2000; do for [i=0:a] {splot 'hermite51-bh' u 2:3:4:1 index i w p ps 1 pt 7 lt palette title (column(1)*100) ; pause 0.01}
def output(t,r,v):
  torb=to()
  for b in range(n):
    print t,
    for c in range(dim):
      print r[b][c],
    for c in range(dim):
      print v[b][c],
    print torb
  print
  print

## Initial Conditions
## Cold-start sphere for a globular cluster of n stars, uniformly distributed
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
        m[cs] *= (1. + np.random.uniform(-mdev,mdev))
      r[cs] = (x, y, z)
      v[cs] = (0, 0, 0)  # Cold Start!
      cs += 1

# times
t=0.     # initial time
tf=20.    # final time (adimensional). Here [G]=[M]=[r]=1 -> [t]=1
dt=0.    # delta t, will be adapted at each step
ts=0.    # auxiliary variable to get the minimum time at each step
tg=0.005   # delta t base, dt = tg * ts 
tout=0.  # when t ~ tout -> print stars position
tsteps=20 # we have an output on intervals tg/tsteps

# Constants
#G=6.67384e-11
G=1.
sigma=1.
dim=3  # 3 space dimensions
R=0.5   # adimensional distance, see note for time
n=3   # number of stars
step=0 # step counter

# Holding vectors
m=np.zeros(n)          # mass
# current
r=np.zeros((n,dim))    # position
v=np.zeros((n,dim))    # velocity
a=np.zeros((n,dim))    # acceleration
k=np.zeros((n,dim))    # jerk, da/dt=d3r/dt3
# previous time (needed for Hermite schema)
rp=np.zeros((n,dim))    # position
vp=np.zeros((n,dim))    # velocity
ap=np.zeros((n,dim))    # acceleration
kp=np.zeros((n,dim))    # jerk, da/dt=d3r/dt3

# Initial conditions
fixmass=True 
mbase=1.
mdev=0.2
initialize()

# change mass for star -> planet earth like, m=1e-6
m[0] = 1.
m[1] = 1.
m[2] = 1e-6
r[0] = (-0.1,0,0)
r[1] = (0.1,0,0)
r[2] = (0.5,0,0)
visviva = np.sqrt(1./la.norm(r))
v[0] = (0, visviva,0)
v[1] = (0,-visviva,0)
visviva = np.sqrt(2./la.norm(r[2]))
v[2] = (0,-visviva,0)

# saving initial energy
e_i=e()

# Algorithm
# Here we go. We need to introduce some changes from 
# leapfrog algorithm, mainly because we have to calculate
# not only a but also k=d3r/dt3
# This introduce an enhancement up to order 4th when compared
# with leapfrog (2nd orden)

# start by calculating a and k
ts = get_ak()
ts = 1.
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
  ts = 1.
  correct() 
  t += dt
  step += 1

# Good bye
# Print final configuration at t=tf
output(t,r,v)
e_f=e()
print "# Initial Total Energy:", e_i
print "# Final Total Energy:", e_f
print "# Relative variation:", (e_f-e_i)/e_i

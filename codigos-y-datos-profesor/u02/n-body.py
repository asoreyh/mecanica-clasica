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

## Clear acceleration
def clear_a():
  for i in range(n):
    a[i]=(0.,0.,0.)

## Gravity 
def update_a():
  for i in range(n):
    for j in range(i+1,n):
      rji=r[j]-r[i]
      aji = G * m[j] * rji / la.norm (rji)**3
      a[i] += aji
      a[j] -= aji

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

# Constants
#G=6.67384e-11
G=1.
dim=3
n=3

# Holding vectors
m=np.zeros(n)          # mass
r=np.zeros((n,dim))    # position
v=np.zeros((n,dim))    # velocity
a=np.zeros((n,dim))    # acceleration

# times
t=0.
tf=1000.
dt=0.01

# Initial conditions

m[0]=2
m[1]=0.002
m[2]=0

r[0] = (0., 0., 0.)
r[1] = (7.5,0,0)
r[2] = (4,0,0)

v[0]=(0,0,0)
v[1]=(0,-1,0)
v[2]=(0,0,0)

e_i=e()

# Algorithm

update_a()

while (t<=tf):
  print t,
  for b in range(n):
    for c in range(dim):
      print r[b][c],
  print

  for i in range(n):
    v[i] += a[i] * dt / 2.
    r[i] += v[i] * dt
  clear_a()
  update_a()
  for i in range(n):
    v[i] += a[i] * dt / 2.
  t += dt

# Good bye
e_f=e()
print "# Energía inicial:", e_i
print "# Energía final:", e_f
print "# Var relativa:", (e_f-e_i)/e_i

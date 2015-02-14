#!/usr/bin/python
# -*- coding: utf8 -*-

import numpy as np

def f(y):
  return r*y*(1.-y)

def g(y):
  return r-2*r*y

r=0.
y=0.2
fv=open("mapeo_{0:04d}.lya".format(int(y*1000)),"w")

for r in np.arange(2.9,4.,0.001):
  n=500
  y=0.2  # cond inicial
  i=0
  while i<300:
    y = f(y)
    i += 1
  i=0
  n=10000
  l=0.
  while i<n:
    y = f(y)
    l += np.log(abs(g(y)))
    i += 1
  l /= n
  fv.write("{0:.5f} {1:.9f}\n".format(r,l))
  print r
fv.close()

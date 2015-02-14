#!/usr/bin/python
# -*- coding: utf8 -*-

import numpy as np

def f(y):
#  return r*y*(1.-y)
  return (k1*y)/((1.+k2*y)**k3)

n=200 # num de ciclos
k1=0.
lim=150
y=2.  # cond inicial
k2=1.
k3=6.
k1=0


for k3 in np.arange(1,7,1):
  fv=open("hassel_{0:04d}.dat".format(int(k3*1)),"w")
  for k1 in np.arange(1.,1000.,0.5):
    y=2.  # cond inicial
    for i in range(n):
      y = f(y)
      if (i>=lim):
        fv.write("{0:.5f} {1:.9f}\n".format(k1,y))
    print k3,k1,y
  fv.close()

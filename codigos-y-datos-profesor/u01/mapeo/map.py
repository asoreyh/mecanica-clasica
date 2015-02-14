#!/usr/bin/python
# -*- coding: utf8 -*-

import numpy as np

def f(y):
  return r*y*(1.-y)

n=100 # num de ciclos
r=3.5475
y=0.2  # cond inicial
fv=open("mapeo_{0:04d}_{1:03d}.dat".format(int(r*1000),int(y*1000)),"w")
for i in range(n):
  y = f(y)
  fv.write("{0:d} {1:.9f}\n".format(i,y))
fv.close()

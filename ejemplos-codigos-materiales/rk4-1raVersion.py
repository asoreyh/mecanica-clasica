#!/usr/bin/python
# -*- coding: utf8 -*-

import numpy as np


def f(t, y):
  return -2.*y+t+4;

def rk4(dt,t,y):
  k1 = dt * f(t, y)
  k2 = dt * f(t + dt/2., y + k1/2.)
  k3 = dt * f(t + dt/2., y + k2/2.)
  k4 = dt * f(t + dt, y + k3)
  return (k1/6.+k2/3.+k3/3.+k4/6.)

dt=0.1
y=1.
t=0.
tf=10.
n=int((tf-t)/dt+0.5)+1

for i in range(n):
  print i,t,y
  y = y + rk4(dt,t,y)
  t += dt

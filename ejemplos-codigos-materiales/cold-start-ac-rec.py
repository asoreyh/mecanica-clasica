#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
cold-start algorithm for the n-body problem
Put n stars uniformily distributed within a 3D-sphere of radius R=1 in
repose.
This version uses acceptance-rejection method. Easier but need to calculate
almost double of points (2^3/(4/3pi) ~ 0.5)
Created on Sun Oct 26 18:51:52 2014
@author: asoreyh
"""

# Numpy!
import numpy as np

R=1.
n=1000

# si fixmass es verdadera, fija las n masas a mbase
# si no, usa una distribución uniforme en mbase (1 +/- mdev)
fixmass=True  # si fixmass, todas las masas son iguales a m
m=1.
mbase=1.
mdev=0.2

# si coldstart es verdadera, todas las velocidades son nulas
# si no, pone una gaussiana en cada componente cartesiana 
# centrada en 0. y con sigma vdev. En este caso hay que restar
# Vcm de todas las velocidades #FIXME
coldstart=True
v=0.
vdev=0.2

i=0
r2=R**2
while (i<n):
  x = np.random.uniform(-1., 1.)
  y = np.random.uniform(-1., 1.)
  z = np.random.uniform(-1., 1.)
  if (x*x+y*y+z*z<=r2):      
    if not fixmass:
      m = mbase * (1. + np.random.uniform(-mdev,mdev))
    
    # coldstart, todas las v son 0
    # TODO, incluir cálculo de velocidades
    if coldstart:
      vx=v
      vy=v
      vz=v
    i += 1
    print i,m,x,y,z,vx,vy,vz

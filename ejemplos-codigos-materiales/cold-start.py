#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
cold-start algorithm for the n-body problem
Put n stars uniformily distributed within a 3D-sphere of radius R=1 in
repose.
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

for i in range(n):
  phi = np.random.uniform(0., 2*np.pi)
  costheta = np.random.uniform(-1., 1.)
  u = np.random.uniform(0., 1.)
  if not fixmass:
    m = mbase * (1. + np.random.uniform(-mdev,mdev))
  
  theta = np.arccos(costheta)
  # esto puede producir problemas si u<0. No es el caso
  r = R * np.power(u,1./3.)  
  
  x = r * np.sin(theta) * np.cos(phi)
  y = r * np.sin(theta) * np.sin(phi)
  z = r * np.cos(theta)
  # coldstart, todas las v son 0
  # TODO, incluir cálculo de velocidades
  if coldstart:
    vx=v
    vy=v
    vz=v

  print i,m,x,y,z,vx,vy,vz

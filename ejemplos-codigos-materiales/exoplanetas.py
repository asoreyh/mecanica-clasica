#!/usr/bin/python
# -*- coding: utf8 -*-

"""
Código hecho en clases prácticas de Introducción a la Física
E. de Física - UIS - Santander - Agosto de 2014

Lee archivo de exoplanetas, Ejemplo sobre como utilizar la nuestra clase vector y el algoritmo visto en clase
para calcular trayectorias en R3
"""

import math
import random

####################### mi clase vector

class vector:
    def __init__(self,lista):
        self.componentes=[]
        for componente_i in lista:
          self.componentes.append(float(componente_i))
        self.dim=len(self.componentes)
        self.mod=math.sqrt(prod_escalar(self,self))

###################### Operaciones entre vectores

def prod_escalar(a,b):
    if (a.dim == b.dim):
        pesc=0.
        for i in range(0,a.dim):
            pesc += (a.componentes[i] * b.componentes[i])
        return pesc
    else:
        print "El productor escalar se define para vectores de la misma dimension"
        return None

def coseno(a,b):
    if (a.dim == b.dim):
      mods=a.mod*b.mod
      if (mods != 0 ):
        return (prod_escalar(a,b) / mods)
      else:
        print "Para calcular el coseno los vectores no pueden ser nulos"
        return None
    else:
        print "Sólo se calcula el coseno para vectores de la misma dimension"
        return None

def suma(a, b):
    suma=[]
    if (a.dim == b.dim):
       for i in range(0,a.dim):
          suma.append(a.componentes[i]+b.componentes[i])
       return vector(suma)
    else:
       print "La suma se define para vectores de la misma dimension"
       return None

def resta(a, b):
    resta=[]
    if (a.dim == b.dim):
       for i in range(0,a.dim):
          resta.append(a.componentes[i]-b.componentes[i])
       return vector(resta)
    else:
       print "La resta se define para vectores de la misma dimension"
       return None

################################## Producto por un escalar

def vector_escalar(a, num):
    escalar=[]
    for i in range(0,a.dim):
      escalar.append(num*a.componentes[i])
    return vector(escalar)

################ Unidades
## Hay que pasar todo a unidades del SI
## Usamos múltiplos comunes (1km = 1000m; 1 dia=86400 s; etc)

dia=86400.
masa_jupiter= 1.898e27
masa_sol= 1.989e30
km=1000.
ua=1.5e8 * km
radio_sol = 695500. * km
radio_jup = 69911 * km
G = 6.67e-11


################ clase exoplaneta
## ejemplo de clase para contener los datos de los planetas que leemos del archivo del catálogo
## haremos una lista de "exoplanetas"
## Si el parámetro no está en el catálogo, pone un valor raro para filtrarlo despues

class exoplaneta:
  def __init__(self,lista):
    if (lista[0] != ''):
        self.nombre=lista[0]
    else:
        self.nombre=None

    if (lista[1] != ''):
        self.semieje=float(lista[1])*ua
    else:
        self.semieje=-1.

    if (lista[2] != ''):
        self.periodo=float(lista[2])*dia
    else:
        self.periodo=-1.
    
    if (lista[3] != ''):
        self.exc=float(lista[3])
    else:
        self.exc=-1.
    
    if (lista[4] != ''):
        self.masa=float(lista[4])*masa_jupiter
    else:
        self.masa=-1.
    
    if (lista[5] != ''):
        self.masae=float(lista[5])*masa_sol
    else:
        self.masae=-1.
    
    if (lista[6] != ''):
        self.radioe=float(lista[6])*radio_sol
    else:
        self.radioe=-1.
    
    if (lista[7] != ''):
        self.tempe=float(lista[7])
    else:
        self.tempe=-1.
    
    if (lista[8] != '\n'):
        self.radio=float(lista[8])*radio_jup
    else:
        self.radio=-1.

################################## mi programa
# aqui comienza mi código...

# Este código tiene dos partes: 
# Primero: lectura y selección de exoplanetas del catálogo
# Segundo: cálculo de la órbita del exoplaneta elegido

# Primera parte

exoplanetas = []
for linea in open("exoplanets.csv"):  # itero sobre todo el archivo problema.dat
  if ((linea.strip()).startswith("#")): # salto los comentarios que empiezan con '#'
    continue
  linea.rstrip() # remueve el salto de línea final
  exo=exoplaneta(linea.split(',')) # creo el objeto planeta de la clase 'exoplaneta' a partir de la lista de parámetro
  #filtro los parámetros que busco: a>0.5 ua, excentricidad > 0.3, y masa de estrella y semieje mayor conocidos:
  if (exo.semieje > 0.5*ua and exo.exc > 0.3 and exo.masae > 0. and exo.masa > 0. and exo.masae/exo.masa > 100.):
    print "#",exo.nombre
    exoplanetas.append(exo) # y si cumple los parámetros de búsqueda lo pongo en la lista de exoplanetas

# ahora exoplanetas es una lista de exoplanetas que cumplen los criterios de selección
# y me quedo con uno al azar: random.randrange(n) devuelve un numero aleatorio entre 0 y n-1, entonces 
# saco un numero aleatorio entre 0 y la cantidad de elementos que tiene la lista de exoplanetas, y uso ese valor
# como índice

myexo = exoplanetas[random.randrange(len(exoplanetas))]

# ahora, myexo es el exoplaneta elegido al azar que cumple con los parámetros de interés
print "# Sobre un total de ",len(exoplanetas)," que cumplen con el criterio de búsqueda"
print "# Se eligió al azar al exoplaneta",myexo.nombre

########## Segunda parte, órbita
# Calculo algunos parámetros orbitales. Trabajo en el SI

foco = myexo.semieje * myexo.exc
afelio = myexo.semieje + foco   # era lo mismo hacer afelio=myexo.semieje*(1.+myexo.exc)
visviva = math.sqrt(G * myexo.masae * ((2./afelio)-(1./myexo.semieje)))
periodo = 2*math.pi*math.sqrt(myexo.semieje**3 / (G*myexo.masae) ) # uso la 3ra ley de kepler para el periodo
# tambien podría haber usado el dato del catálogo: myexo.periodo
intervalos = 1000 # en 1000 pasos completo una órbita
dt=periodo / intervalos 
n = 1 # numero de orbitas a recorrer

##### CONDICIONES INICIALES
# Trabajo en 2D porque por Kepler las órbitas están en un plano
# Necesito el vector posición del planeta. El origen está en el foco donde está la estrella.
# posicion inicial, empiezo en el afelio
r=vector([-afelio,0])
# ahora el vector velocidad, perpendicular a r, lo pongo en la dirección y y uso visviva
v=vector([0,visviva])  # si uso -visviva cambio el sentido de giro
# y defino el vector aceleración, que por ahora vale 0 (se calcula en el bucle)
a=vector([0,0]) # empiezo con 0, porque despues la cambio 

# Dos auxiliares
r_unit=vector([0,0])  # lo uso para el vector unitario
amod = 0. # será el módulo de la aceleración

for i in range(0, int(n*intervalos)):
    # imprimo el paso i, el tiempo transcurrido en días y la posición del planeta (en unidades astronomicas) 
    print i,i*dt/dia ,r.componentes[0]/ua ,r.componentes[1]/ua

    # Algoritmo Newton-Hooke
    # actualizo la posición r_{i+1} = r_i + dt * v_i 
    r=suma(r,vector_escalar(v,dt))
    # determino un vector con dirección r, sentido hacia dentro y módulo 1
    if (r.mod):
      r_unit=vector_escalar(r,(-1.)/r.mod)
    else:
      print "Hay un problema: |r| = 0"
      exit()

    # ahora, calculo la aceleración
    # en kepler, la aceleración centrípeta está dada por la gravedad 
    # que tiene dirección radial, sentido hacia dentro y módulo 
    # |a| = G * M / |r|^2
    amod = G * myexo.masae / (r.mod**2)
    # entonces la aceleración es:
    # a = amod * r_unit 
    a=vector_escalar(r_unit,amod)
    # y actualizo la velocidad v_{i+1} = v_i + a_i dt
    v=suma(v,vector_escalar(a,dt))
    # e itero 

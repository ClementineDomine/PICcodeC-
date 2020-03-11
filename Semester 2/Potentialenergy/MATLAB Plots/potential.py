# -*- coding: utf-8 -*-
# Written by Clementine Domine
# Last modified on 22/02/0202
# Plot a Phase space animation of an arbitrary number of particles in a 1d plasma at fixed radius

import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation
import math

#----------------------------------------------
#INPUT:
row1= [];
row2= [];
row3= [];
row4= [];
row5= [];

phic1='potentialenergycorrection.csv'
phist1='potentialenergyselftot.csv'
phisc1='potentialenergyselfcorrect.csv'

phitc1='potentialenergytotalcorrect.csv'
phitrap1='potentialenergytrap.csv'
#----------------------------------------------
phic = list(csv.reader(open(phic1), delimiter=','))
for row in phic:
  for i in range(len(row)):
        row[i] = float(row[i])
        row1.append(row[i])
        
phist = list(csv.reader(open(phist1), delimiter=','))
for row in phist:
  for i in range(len(row)):
        row[i] = float(row[i])
        row2.append(row[i])
        
phisc = list(csv.reader(open(phisc1), delimiter=','))
for row in phisc:
  for i in range(len(row)):
        row[i] = float(row[i])
        row3.append(row[i])
    
phitc = list(csv.reader(open(phitc1), delimiter=','))
for row in phitc:
  for i in range(len(row)):
        row[i] = float(row[i])
        row4.append(row[i])



phitrap = list(csv.reader(open(phitrap1), delimiter=','))
for row in phitrap:
    for i in range(len(row)):
        row[i] = float(row[i])
        row5.append(row[i])
length = [i for i in range(0,10000)] 
plt.figure(num=None, figsize=(9, 7), dpi=80, facecolor='w', edgecolor='k')

plt.scatter(length,row5,s=0.1,label='U self correction ')
plt.scatter(length,row1,s=0.1,label='U self total')
plt.scatter(length,row3,s=0.1,label='U self correct')
plt.scatter(length,row4,s=0.1,label='UW +US correct')
plt.scatter(length,row2,s=0.1,label='U Well total')



plt.ylim(-0.7*10 **(-15), 0.1*10 **(-15))
#ax = matplotlib.pyplot.axes(xlim = , ylim = )
plt.title("Potential energy per particles")
plt.legend()
plt.ylabel("Potential energy (V/particles)")
plt.xlabel("Z (m)")
#matplotlib.pyplot.tight_layout()
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0), useOffset=False, useLocale=False, useMathText=True)
plt.show()
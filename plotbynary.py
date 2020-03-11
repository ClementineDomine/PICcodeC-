#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  7 16:08:48 2020

@author: clementine
"""

import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation
import math
from astropy.io import ascii

from sys import argv



f=open("herebin.bin","rb")

x = np.fromfile(f)
print (x)

#num=list(f.read())
#print (num)
f.close()


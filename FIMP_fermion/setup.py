import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.special import kn
from scipy.interpolate import interp1d
import astropy.constants as ct
import pandas as pd
import subprocess
import sys
import os

from FIMP_fermion import *

vals = np.logspace(-3,3)
model = FIMP_fermion(vals,2000.,125.,0.1,1e-9,1e-9)
#model.plot_omega(False,True,-3,2.99,100,-12,-9,4,colors=['r','g','b','darkviolet','pink','chartreuse','deepskyblue','plum'])
model.plot_y(-3,4,100,1,4,4,['r','g','b','darkviolet'])
#print(model.relic_abundance(small_coupling=False,micromegas=True))
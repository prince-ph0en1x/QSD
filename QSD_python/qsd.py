## \date: 14-03-2018 -
## \repo: https://gitlab.com/prince-ph0en1x/QSD
## \proj: OpenQL Libraries for Q Algo. development

###################################################################################################
#                                            Algorithm                                            #
###################################################################################################
    
# Quantum Shannon Decomposition

# Reference: Synthesis of Quantum Circuits (Shende et.al)

###################################################################################################

from openql import openql as ql
from qxelarator import qxelarator
from random import random
#from math import *
from cmath import *
import numpy as np
import os

###################################################################################################

# Notes/Snippets

'''
cd /mnt/7A06EEA206EE5E9F/GoogleDrive/TUD_CE/Thesis/SimQG/QSD/QSD_python
python3 qsd.py
../QXSim/qx_simulator_1.0.beta_linux_x86_64 test_output/qg.qasm
'''

'''
qx = qxelarator.QX()
qx.set('basic.qasm')
qx.execute()
res = qx.get_measurement_outcome(0)
'''

###################################################################################################

im = complex(0.0,1.0)

###################################################################################################

# Quantum Gates

def AP(t):
    return np.matrix([[exp(im*t), 0], [0, exp(im*t)]])

def Rx(t):
    return np.matrix([[cos(t/2), im*sin(t/2)],[im*sin(t/2), cos(t/2)]])

def Ry(t):
    return np.matrix([[cos(t/2), sin(t/2)],[-sin(t/2), cos(t/2)]])

def Rz(t):
    return np.matrix([[exp(im*t/2), 0], [0, exp(-im*t/2)]]) 
    # TBD: Check conventional sign of matrix

H = (AP(-pi/2) @ Ry(pi/2) @ Rx(pi)).real
X = (AP(-pi/2) @ Rx(pi)).real
I = np.matrix([[1, 0], [0, 1]]) 
O = np.matrix([[0, 0], [0, 0]])
CX = (np.concatenate((np.concatenate((I,O)) ,np.concatenate((O,X)) ), axis=1)).real
XC = np.kron(H,H) @ CX @ np.kron(H,H)
SWAP = CX @ XC @ CX

###################################################################################################


print(XC)
#!/usr/bin/env python
import numpy as np
from scipy.linalg import expm

A = np.array([[0.89495215, -0.80105662, 0.51002593, 0.82839079, 0.34766935],
              [-0.58118727, 0.21694723, -0.22004884, -0.053947705, -0.70431062], 
              [-0.87149079, 0.060792935, -0.52037241,  0.3570078, -0.28580946], 
              [-0.25978119, -0.034334584, -0.6887698, 0.60344181, 0.23295322], 
              [0.97005764, 0.65041857, 0.44865746, -0.083630837, -0.11497899]] )
print(A)
print(expm(1.23*A))
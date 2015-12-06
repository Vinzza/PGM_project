#! /usr/bin/env python
# -*- encoding : utf-8 -*-

import numpy as np
import scipy.linalg as ln
import numpy.random as rd

def whit(x):
	C = np.cov(x)
	w,V = ln.eigh(C)
	A = V.dot(np.diag(w**(-0.5)).dot(V.T))
	x = A.dot(x)
	return x

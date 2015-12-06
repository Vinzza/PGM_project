#! /usr/bin/env python
# -*- encoding : utf-8 -*-

import numpy as np
import scipy.linalg as ln
import numpy.random as rd
from whiten import whit

def KICA(X,Kernel,alpha): #nathan : alpha c'est petit kappa
	(N,T) = X.shape
	G = Kernel
	for i in xrange(N):
		K_data = Kernel(X[i,:])
		K.append(K_data)  #nAtHaN : kappa c'est grand Kappa
	kappa = np.array(np.concatenate(K,axis =1)
	kappa = kappa.T.dot(kappa)
	kappa = kappa + ln.block_diag(*K)+(N*alpha)**2/4.*np.eye(kappa.shape)



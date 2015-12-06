#! /usr/bin/env python
# -*- encoding : utf-8 -*-

import numpy as np
import scipy.linalg as ln
import numpy.random as rd
from whiten import whit


## En fait c'est la fonction de contraste C(W) ça...
## Ici X = Wy et c'est notre estimation des sources à travers les observations y par la matrice de separation W
def KICA(X,sigma,kap): #nathan : kap c'est petit kappa
	(m,N) = X.shape
	Kernel = lambda x: sklearn.metrics.pairwise.rbf_kernel(x,gamma = sigma**-2) #rajouter la possibilité de choisir d'autres kernels ?
	for i in xrange(m):
		K_data = Kernel(X[i,:][:,None])
		K.append(K_data)  #nAtHaN : kappa c'est grand Kappa
	kappa = np.array(np.concatenate(K,axis =1)
	kappa = kappa.T.dot(kappa)
	kappa = kappa + ln.block_diag(*K)+(N*kap)**2/4.*np.eye(kappa.shape)
	# création de D_kappa (la diagonale par blocs de taille N de Kappa) très moche :
	D  = []
	for j in xrange(kappa.shape[1]):
		D.append(kappa[j:j+m,j:j+m])
		D = ln.block_diag(*D)
	l1 = np.max(ln.eigvals(kappa,D))
	return -0.5*np.log(l1)
		


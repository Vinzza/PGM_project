#! /usr/bin/env python
# -*- encoding : utf-8 -*-

import numpy as np
import time
import scipy.linalg as ln
import numpy.random as rd
import sklearn.metrics.pairwise as skl


## En fait c'est la fonction de contraste C(W) ca...
## Ici X = Wy et c'est notre estimation des sources a travers les observations y par la matrice de separation W
def contrast(X,sigma,kap): #nathan : kap c'est petit kappa
	(m,N) = X.shape
	Kernel = lambda x: skl.rbf_kernel(x,gamma = sigma**-2) #rajouter la possibilite de choisir d'autres kernels ?
	K = []
	for i in xrange(m):
		K_data = Kernel(X[i,:][:,None])
		K.append(K_data)  #nAtHaN : kappa c'est grand Kappa
	a = time.time()
	kappa = np.array(np.concatenate(K,axis =1))
	kappa = kappa.T.dot(kappa)
	kappa = kappa + ln.block_diag(*K)+(N*kap)**2/4.*np.eye(*kappa.shape)
	b = time.time()
	print b-a
	# creation de D_kappa (la diagonale par blocs de taille N de Kappa) tres moche :
	D  = []
	for j in xrange(m):
		D.append(kappa[j*N:(j+1)*N,j*N:(j+1)*N])
	D = ln.block_diag(*D)
	a = time.time()
	l1 = np.max(np.real(ln.eigvals(kappa,D)))
	b = time.time()
	print b-a
	return -0.5*np.log(l1)

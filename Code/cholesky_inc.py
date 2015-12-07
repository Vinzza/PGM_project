#! /usr/bin/env python
# -*- encoding : utf-8 -*-

import numpy as np
import scipy.linalg as ln

def ichol(K,eta):
	N = K.shape[0]
	if N != K.shape[1]:
		print "La matrice d'entree n'est pas carree"
		return
	i = 0
	P = np.eye(N,N)
	G = np.diag(np.diag(K))
	K1 = K.copy()
	while np.sum(np.diag(G)) > eta and i < N:
		j = np.argmax(np.diag(G)[i:])
		P[j,j], P[i,i],	P[i,j], P[j,i] = 0.,0.,1.,1.
		K1[:,i],K1[:,j] = K1[:,j],K1[:,i]
		K1[i,:],K1[j,:] = K1[j,:],K1[i,:]
		G[i,:i+1],G[j,:i+1] = G[j,:i+1],G[i,:i+1]
		G[i,i] = np.sqrt(K1[i,i])
		G[i+1:,i] = (K1[i+1:,i]- np.sum([G[i+1:,j]*G[i,j] for j in xrange(i)],axis=0))/G[i,i]
		for k in xrange(i+1,N):
			g = K[k,k] - np.sum(G[k,:i+1]**2)
			G[k,k] = g
		print np.sum(np.diag(G))
		i += 1
	return P,G,i

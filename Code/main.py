#! /usr/bin/env python
# -*- encoding:utf-8 -*-

from whitening import whiten
import sklearn as skl


## Data creation

## Data treatment
X = X -np.mean(X,axis=1)[:,None]
X = whiten(X)

## Kernel definition



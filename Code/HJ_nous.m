function [G,x] = HJ(y)
%% G is the mixing matrix, x the sources, and y the observations, such that we have y =G.x
%% y is a matrix of n column (sources) and T lines for the temporality

% Data conditionning
[T,n] = size(y);
y = y - mean(y);

W = rand(n,n);
Wold = rand(n,n);
eps = 10^-3;
f = lambda x: 
g = lambda x: 

while norm(W-Wold) > eps


end
end

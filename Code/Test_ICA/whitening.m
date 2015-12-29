function [X_white, mu, white_matrix] = whitening( X )

m = size( X, 1 );
N = size( X, 2 );

epsilon = 0.0001;

% removing the mean of all column
mu = mean(X,2);
X = bsxfun(@minus, X, mu);

% Calculating the Covariance Matrix
Cov = X * X';
[V,D,~] = svd(Cov);

% Whitening of the data
D_dum = diag(sqrt( 1 ./ (diag(D)+epsilon) ));
white_matrix = V * D_dum * V';
X_white = white_matrix * X;

end
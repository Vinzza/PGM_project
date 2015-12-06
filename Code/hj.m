function [B, mdS] = hj(X, k, varargin)
% HJ    Use the Herault-Jutten (HJ) algorithm to learn the weights of the
% matrix b for ICA.
%
% B = hj(X, options)
%   X is the set of input vectors.  Each column of X is one sample.
%   k is the number of independent components to extract.
%   B is the seperating matrix, Y = BX are the independent signals
%
% options (specified by key/value pairs)
%   'rate' = 0.5        the learning rate (initial)
%   'B' = rand          the initial seperating matrix
%   'niter' = 1         the number of iterations to run
%   'g' = t.^3          the 'g' independence forcing function
%   'h' = t             the 'h' independence forcing function

%
% David Gleich
% CS 152 - Neural Networks
% 12 December 2003
%

% dimensionality assessment
[m n] = size(X);
if (k > m)
    error('Cannot extra more sources than sensors.');
end;

% subtract the mean
X = X - repmat(mean(X,2), 1, n);

% options
options = struct(...
    'rate', 0.5, ...
    'S', rand(k,m), ...
    'niter', 1, ...
    'g', inline('t.^3'), ...
    'h', inline('t'));
options = getopt(options, varargin);

S = options.S;
u = options.rate;
g = options.g;
h = options.h;

% make sure the S matrix has zero diagonal
S = S - diag(diag(S));

% compute B "exactly"
%B = eye(k,m) - S;
B = (eye + S)^-1;

success = 0;
mdS = 0;

for niter=1:options.niter
    Bold = B;
    % update B, S
    for ii = 1:n
        y = B*X(:, ii);
        dS = u*g(y)*h(y');
        mdS = [mdS norm(dS, 'fro')];
        S = S + dS;
        S = S - diag(diag(S));
        B = (eye + S)^-1;
    end;
    if (success == 1)
        break;
    end;
end;

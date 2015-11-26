function [B] = easi(X, k, varargin)
% EASI    Use the EASI algorithm to learn the weights of the
% matrix B for ICA.
%
% B = easi(X, k, options)
%   X is the set of input vectors.  Each column of X is one sample.
%   k is the number of independent components to extract.
%   B is the seperating matrix, Y = BX are the independent signals
%
% options (specified by key/value pairs)
%   'rate' = 0.001      the learning rate (initial)
%   'B' = rand          the initial seperating matrix
%   'niter' = 1         the number of iterations to run
%   'g' = t             the 'g' independence forcing function
%   'h' = tanh(t)       the 'h' independence forcing function

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
    'rate', 0.001, ...
    'B', rand(k,m), ...
    'niter', 1, ...
    'g', inline('t'), ...
    'h', inline('tanh(t)'));
options = getopt(options, varargin);

B = options.B;
u = options.rate;
g = options.g;
h = options.h;

success = 0;
iter = 1;

for niter=1:options.niter
    Bold = B;
    
    % update B, S
    for ii = 1:n
        y = B*X(:, ii);
        B = B - u*(y*y' - eye(k,k) + g(y)*h(y') - h(y)*g(y'))*B;
        
        if (any(~isfinite(B)))
            warning(sprintf('Lost convergence at iterator %i; lower learning rate?', iter));
            success = 11;
            break;
        end;
        
        iter = iter + 1;
    end;
   
    
    if (success == 1)
        break;
    end;
end;
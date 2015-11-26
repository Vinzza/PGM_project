function w=oja(X, varargin)
% OJA  Oja's rule for determining the weights of a principal component
%   neural network.
%
% w = oja(X, options)
%   X is the set of input vectors.  Each column of X is one sample.  
%   w is the first principal component
%
% options (specified by key/value pairs)
%   'rate' = 0.1        the learning rate
%   'w' = rand          the initial vector
%   'niter' = 1         the number of iterations
%   'delta' = 0.001     the convergence criteria
%   'animate' = 0       display an animation of the sucessive vectors


%
% David Gleich
% CS 152 - Neural Networks
% 12 December 2003
%

% get the dimensionality
[m n] = size(X);

% random initial weights
w = rand(m,1);

options = struct( ...
    'rate', .1, ...
    'w', w, ...
    'niter', 1, ...
    'delta', .001, ...
    'animate', 0);
options = getopt(options, varargin);

% validate options
if (size(options.w,1) ~= m)
    warning('Invalid initial set of weights, using random.');
    options.w = w;
end;

w = options.w;

% check for animation
if (options.animate > 0)
    w = animated_oja(X, options);
    return;
end;
    
success = 0;

% run through all input samples
for iter = 1:options.niter
    wold = w;
	for ii = 1:n 
        y = w'*X(:,ii);
        % y is a scalar, not a vector
        w = w + options.rate*(y*X(:,ii) - y^2*w);
    end;
    
    if (any(~isfinite(w)))
        warning(sprintf('Lost convergence at iterator %i; lower learning rate?', iter));
        success = 11;
        break;
    end;
    
    if (norm(w - wold) < options.delta)
        success = 1;
        break;
    end;
    
    if (success > 0)
        break;
    end;
end;

if (success == 0)
    warning(sprintf('Did not converge in %i iterations: possibly unstable input.', options.niter));
end;

% ANIMATED_OJA An animated implementation of Oja's rule.
%   All non-standard options must be set in the base oja function.  
function w = animated_oja(X, options)

% get the dimensionality
[m n] = size(X);

if (m ~= 2)
    error('Animation only works for two dimensional signals.');
end;
    
w = options.w;
success = 0;

% run through all input samples
for iter = 1:options.niter
   
    wold = w;
	for ii = 1:n 
        
        % plot current results
        plot(X(1,:), X(2,:), '.', [0 w(1)], [0 w(2)], '-', 'LineWidth', 3);
        drawnow;
        
        y = w'*X(:,ii);
        % y is a scalar, not a vector
        w = w + options.rate*(y*X(:,ii) - y^2*w);
    end;
    
    if (norm(w - wold) < options.delta)
        % plot last result
        success = 1;
        break;
    end;
    
    if (success == 1)
        break;
    end;
end;

if (success ~= 1)
    warning(sprintf('Did not converge in %i iterations: possibly unstable input.', options.niter));
end;

plot(X(1,:), X(2,:), '.', [0 w(1)], [0 w(2)], '-', 'LineWidth', 2);

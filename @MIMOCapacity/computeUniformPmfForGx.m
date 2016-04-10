function [pmf, reachable] = computeUniformPmfForGx(G, x_max, ...
    y_min, delta, nbins)
% computeUniformPmfForGx computes the PMF for a uniformly distributed
% random variable y, where y = G*x where x is bounded to the interval [0,
% x_max].  
%
% We're using the channel matrix y = G*x + w.
%
%   pmf = computeUniformPmfForGx(G, x_max, y_min, delta, nbins)
%
% pmf is the probability mass function for a random variable that is
%   uniformly distributed over the range of values G*x for for x(d) in the
%   interval [0, x_max(d)] for all dimensions d.  pmf is a is a
%   multidimension matrix of size nbins.  
% reachable is a multidimension matrix of the same size as pmf that
%   indicates whether the bin of the pmf is reachable by G*x for the
%   bounded values of x.  This output is provided for debugging and testing
%   purposes.  
%
% G (n_r*n_t matrix) is the channel matrix for n_t transmitters and n_r
%   receivers.  
% x_max (scalar or n_t-element column vector) is the maximum signal that
%   can be transmitted by each transmitter. 
% y_min (n_r-element row vector) is the lowest y (or G*x) that is included
%   in a bin in the PMF.
% delta (n_r-element row vector) is the size (in terms of y) of each bin of
%   the PMF. nbins (n_r-element row vector) is the number of bins along
%   each dimension of y (or of the PMF) in the PMF. 

%% check inputs and get dimensions

% G
if(~ismatrix(G))
    error('Parameter G must be 2D matrix.');
end
[n_r, n_t] = size(G);

% x_max
if(isscalar(x_max))
    x_max = repmat(x_max, n_t, 1);
end
if(~isequal(size(x_max), [n_t, 1]))
    error(['Parameter x_max must be a n_t-element column vector ' ...
        '(or a scalar).']);
end
if(any(x_max<0))
    error('Parameter x_max should be non-negative.');
end

% y_min
if(~isequal(size(y_min), [1, n_r]))
    error('Parameter y_min must be a n_r-element row vector.');
end
if(any(y_min>0))
    warning('The PMF excludes the y = 0 point.');
end

% delta
if(~isequal(size(delta), [1, n_r]))
    error('Parameter delta must be a n_r-element row vector.');
end
if(any(delta <= 0))
    % delta should not be equal zero; otherwise, we'd have duplicate bins
    % covering the same values (assuming that we have more than one bin in
    % that dimension).
    error('Parameter delta must be positive.');
end

% nbins
if(~isequal(size(nbins), [1, n_r]))
    error('Parameter nbins must be a n_r-element row vector.');
end
if(any(nbins < 1))
    error('Parameter nbins should be a whole number.');
end

% check that G*x_max is included in the PMF
Gx_max = G*x_max; % column vector
y_max = y_min + nbins .* delta; % row vector
if(any(Gx_max' > y_max))
    warning('The PMF excludes the y = G*x_max point.');
end

%% determine whether each bin of the PMF is reachable

% indexing weights for matrix pmf for MIMOCapacity.convertToLinearIndex
iw = MIMOCapacity.convertToLinearIndexWeights(nbins);

% Pre-allocate a multi-dimension matrix to store whether a bin of the PMF
% is reachable with G*x given bounded x.
% The last dimension of 1 in the argument of zeros ensures that zeros does
% not return a square matrix when n_r is 1 (and doesn't do anything
% otherwise).
reachable = false([nbins, 1]);

% precompute the parameters of linprog that don't vary between bins
f = zeros(1, n_t);
A = [G; -G ];
lb = zeros(1, n_t);
ub = x_max';

% optimization options
% As recommended by the MATLAB documentation at
% http://www.mathworks.com/help/optim/ug/choosing-a-solver.html, use the
% 'interior-point' instead of the default 'interioir-point-legacy'
% algorithm.
% TODO: After I understand linear programming algorithms more, I may be
% able to choose options to reduce the computation time.
lpopt = optimoptions('linprog', 'Algorithm', 'interior-point', ...
    'Display', 'off');

% fill in each bin of pmf in parallel
parfor ibin = 1:prod(nbins) % ibin is the linear index for the PMF bin
    
    % Get the matrix subscript index so we can determine the bounds of this
    % bin.
    ms = MIMOCapacity.convertLinearToSubscriptIndex(iw, ibin);
    
    % lower bounds and upper bounds for this bin.
    y_lb = y_min + (ms-1) .* delta;
    y_ub = y_min + ms .* delta;
    
    % y_ub and y_lb are row vectors; b should be a column vector.
    b = [y_ub'; -y_lb'];
    
    [~, ~, exitflag] = linprog(f, A, b, [], [], lb, ub, [], lpopt);
    
    switch exitflag
        case 1
            % Function converged to a solution x.
            reachable(ibin) = true;
        case -2
            % No feasible point was found.
            reachable(ibin) = false;
        otherwise
            keyboard
            error('Unexpected exitflag %d from linprog for ibin %d.', ...
                exitflag, ibin);
    end % switch exitflag
    
end % parfor ibin

%% compute the uniform PMF
% Here, we rely on the approximation that all reachable bins can be reached
% with the same probability.  In actuality, some of the bins at the outer
% boundary of the reachable values of G*x may have a lower probability than
% those bins that are entirely witin the reachable range of G*x.  
% The error due to this error should decrease as nbins increases.  

pmf = zeros([nbins, 1]);
pmf(reachable) = 1/nnz(reachable);

end


function pmf = generateUniformPmfForGx(G, x_max, ymin, delta, nbins)
% generateTransformedUniformPmf generates the uniformly-distributed PMF
% that covers the possible values of G*x where x in each dimension is in
% the range [0, x_max].  
%
% Uses the channel model: y = G*x + w.
%
%   PMF = generateTransformedUniformPmf( ...
%       G, X_MAX, YMIN, DELTA, NBINS)
%
% PMF (n_r-dimension matrix) is the generated PMF.
%
% G (n_r*n_t matrix) is the channel matrix.  
% X_MAX (scalar or n_t-element column vector) is the maximum signal that
%   can be transmitted by each transmitter.
% YMIN (n_r-element row vector) is the minimum value (respectively) covered
%   covered by the PMF along each dimension. 
% DELTA (n_r-element row vector) is the size (in terms of y) of each bin of
%   the PMF.
% NBINS (n_r-element row vector) is the number of bins along each dimension
%   of y.  NBINS should consist of whole numbers.  
%
% Note that nbins, ymin, and ymax should be specified to be large enough to
% cover G*x (without the noise from w).  When the PMF of G*x is convolved
% against the PMF of the noise w, the matrix would be correspondingly
% expanded (using the same bin sizes).  

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
    error(['Parameter X_MAX must be a n_t-element column vector ' ...
        '(or a scalar).']);
end
if(any(x_max<0))
    error('Parameter X_MAX should be non-negative.');
end

% ymin
if(~isequal(size(ymin), [1, n_r]))
    error('Parameter YMIN must be a n_r element row vector.');
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
    error('Parameter NBINS must be a n_r element row vector.');
end
if(any(nbins<1))
    error('PARAMETER NBINS should be a whole number.');
end

% check that G*x_max is included in the PMF
Gx_max = G*x_max; % column vector
y_max = y_min + nbins .* delta; % row vector
if(any(Gx_max' > y_max))
    warning('The PMF excludes the y = G*x_max point.');
end

%% Check that G*x_max is at least close to the maximum-value corner of its
% bin.  Otherwise, it may be very improbable to hit the bin that contains
% G*x_max.

% Transpose Gx_max for convertPointToSubscriptIndex, which expects each
% point to be represented by a row in input y.  
index_Gx_max = MIMOCapacity.convertPointToSubscriptIndex(Gx_max', ...
    ymin, delta, nbins);

% subscript index of y=0
index_zero = MIMOCapacity.convertPointToSubscriptIndex(zeros(1, n_r), ...
    y_min, delta, nbins);

% lower bound of the bin that contains Gx_max
lb_bin_Gx_max = y_min + (index_Gx_max-1) .* delta;

% Position of Gx_max within its bin along each dimension, normalized so
% that 0 means Gx_max is at the lower bound of the bin and 1 means that
% Gx_max is at the upper bound of the bin.  
normalized_pos_of_Gx_max_within_bin = (Gx_max' - lb_bin_Gx_max) ./ delta;

% If in any dimension of y, the bin index_Gx_max does not contain the
% entire range of possible G*x, then ensure that Gx_max is near the upper
% bound of the bin.  Otherwise, if the distribution of G*x barely touches
% bin index_Gx_max, it may take a very long time for the Monte Carlo
% simulation to acquire enough hits in bin index_Gx_max.  
if(any(index_Gx_max > index_zero & ...
        normalized_pos_of_Gx_max_within_bin < 0.9))
    warning(['Gx_max does not extend far enough into its bin; ' ...
        'the probability of hitting this bin may be low, which may ' ...
        'cause this Monte Carlo simulation to take a very long time ' ...
        'to complete.']);
end

%% initial setup

% The number of trials to run together through vectorization.
% Select a value that's
%  - big enough to fully/efficiently utilize the CPU(s),
% - is divisible in many ways, and
% - is not too big (to keep memory use small and to avoid extending the
%   computation time too long). 
trials_per_batch = 1474560; % 45*2^15; 11.25 MiB per double vector.

% The threshold minimum count (to tell when we've probably filled in every
% possible value that can be reached).
% lab book #4, p. 42-43.
tc = 64 + log(prod(nbins));

% Allocate matrix to count hits in each of the bins.
% The last dimension of 1 in the argument of zeros ensures that zeros does
% not return a square matrix when n_r is 1 (and doesn't do anything
% otherwise).
hits = zeros([nbins, 1]);

% pre-compute the weights needed to convert a row of index_x into a linear
% index (needed to index elements in arbitrary-dimension matrices).
indexing_weights = MIMOCapacity.convertToLinearIndexWeights(nbins);

% Pre-compute the matrix index for the bin corresponding to G*x for the
% maximum value of x.  This way, it's more efficient to check whether we
% have the requisite tc hits in this bin.  
li_Gx_max = MIMOCapacity.convertToLinearIndex(indexing_weights, ...
    index_Gx_max);

if(isempty(li_Gx_max))
    error(['Unable to calculate the linear index for the PMF bin ' ...
        'that contains Gx_max.']);
end % if(isempty(li_Gx_max))

%% include x_max into G
% rand() gives uniformly distributed random values in the interval (0,1).
% However, we want x to be uniformly distributed in the interval (0, x_max)
% for each dimension.  
% We generate many many values for x; rather than multiply each one by
% x_max to get the appropriate range of random values, multiple each column
% of G by the x_max for that transmitter.  This should yield the same
% result as multiplying each column of x by the corresponding column of
% x_max before multiplying by G.  

% x_max a column vector; transpose it to get a row vector.
G = bsxfun(@times, G, x_max');

%% count hits until we have at least tc in the bin for Gx_max

% Note: waitbar should be removed before running on the cluster (with
% -nodisplay) to avoid warnings.  
wb = waitbar(0, 'Starting...');

while(hits(li_Gx_max) < tc)
    waitbar(hits(li_Gx_max)/tc, wb, ...
        sprintf('Gxmax bin reached %u of %.0f hits', ...
        hits(li_Gx_max), tc));
    
    nestedTallyGx
end % while(hits(li_Gxmax) < tc)

close(wb)

if(nnz(hits) == 0)
    error('MIMOCapacity:generateTransformedUniformPmf:NoHits', ...
        ['No hits were recorded. The specified YMIN and YMAX do ' ...
        'not adequately cover the possible values of G*x.']);
end % if(nnz(hits) == 0)

%% return a uniform PMF that reaches the same values

% Probability for each bin hit.
unif_prob_per_bin = 1/nnz(hits);

% Generate a PMF that is unif_prob_per_bin where G*x reaches and 0
% elsewhere.
pmf = zeros([nbins, 1]);
pmf(hits>0) = unif_prob_per_bin;

%% nested function to run batch of Monte Carlo simulations
% Since we run this same code within two separate while loops, put it into
% a nested function.  

    function nestedTallyGx
        
        % Sample normalized values of x (since we already included the
        % scaling for x_max in G).
        x = rand(trials_per_batch, n_t);

        % Multiply this way to have each row of Gx be a point.
        Gx = x * G';    

        % Convert each of the points in Gx into a corresponding linear
        % index (li_Gx) for matrix hits.
        msi_Gx = MIMOCapacity.convertPointToSubscriptIndex(Gx, ...
            ymin, delta, nbins);
        li_Gx = MIMOCapacity.convertToLinearIndex( ...
            indexing_weights, msi_Gx);

        % Tally each hit in matrix hits
        for ii = 1:length(li_Gx)
            index_bin = li_Gx(ii);
            hits(index_bin) = hits(index_bin) + 1;
        end % for ii = 1:length(li_Gx)
        
    end % function nestedTallyGx

end

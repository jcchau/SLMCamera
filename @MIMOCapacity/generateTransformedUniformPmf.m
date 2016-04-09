function pmf = generateTransformedUniformPmf( ...
    G, x_max, nbins, ymin, ymax)
% generateTransformedUniformPmf generates the uniformly-distributed PMF
% that covers the possible values of G*x where x in each dimension is in
% the range [0, x_max].  
%
% Uses the channel model: y = G*x + w.
%
%   PMF = generateTransformedUniformPmf( ...
%       G, X_MAX, NBINS, YMIN, YMAX)
%
% PMF (n_r-dimension matrix) is the generated PMF.
%
% G (n_r*n_t matrix) is the channel matrix.  
% X_MAX (scalar or n_t-element column vector) is the maximum signal that
%   can be transmitted by each transmitter.
% NBINS (n_r-element row vector) is the number of bins along each dimension
%   of y.  NBINS should consist of whole numbers.  
% YMIN and YMAX (each a n_r-element row vector) are the minimum and
%   maximum values (respectively) covered covered by the nbins bins along
%   each dimension.  
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

% nbins
if(~isequal(size(nbins), [1, n_r]))
    error('Parameter NBINS must be a n_r element row vector.');
end
if(any(nbins<1))
    error('PARAMETER NBINS should be a whole number.');
end

% ymin and ymax
if(~isequal(size(ymin), [1, n_r]))
    error('Parameter YMIN must be a n_r element row vector.');
end
if(~isequal(size(ymax), [1, n_r]))
    error('Parameter YMAX must be a n_r element row vector.');
end
if(any(ymin >= ymax))
    error('Parameter YMAX should be greater than YMIN.');
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

% The size of a bin along each dimension
delta = (ymax - ymin) ./ nbins;

% Pre-compute the matrix index for the bin corresponding to G*x for the
% maximum value of x.  This way, it's more efficient to check whether we
% have the requisite tc hits in this bin.  
Gx_max = G * x_max;
% Transpose Gx_max for convertPointToSubscriptIndex, which expects each
% point to be represented by a row in input y.  
index_Gx_max = MIMOCapacity.convertPointToSubscriptIndex(Gx_max', ...
    ymin, delta, nbins);
li_Gx_max = MIMOCapacity.convertToLinearIndex(indexing_weights, ...
    index_Gx_max);

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

if(isempty(li_Gx_max))
    warning('MIMOCapacity:generateTransformedUniformPmf:GxmaxOutside', ...
        ['The YMIN and YMAX specified does not include G*x where x ' ...
        'is the maximum as specified by X_MAX. ' ...
        'Skipped running trials until the bin for Gx_max has at ' ...
        'least tc hits.']);
    
    % run one batch so we don't skip the next step
    nestedTallyGx
end % if(isempty(li_Gx_max))

%% and until we have at least tc hits in the bin with fewest non-zero hits

while(min(hits(hits>0)) < tc)
    waitbar(min(hits(hits>0))/tc, wb, ...
        sprintf('Minimum hit count reached %u of %.0f hits', ...
        min(hits(hits>0)), tc));
    
    nestedTallyGx
end % while(min(hits(hits>0)) < tc)

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

function [C, n, t] = computeTableOfCapacity(A, nStart)
% computeTableOfCapacity computes a table of capacity as a function of A.
%
%   [C, n, t] = computeTableOfCapacity(A)
%
% C is the capacity as computed via Smith1971.
% n is the optimum number of discrete constellation points in X.
% t is the incremental time needed for the capacity computation.
%
% A (vector) is the amplitude limit of X (in Y=X+N).  A should be sorted in
%   ascending order.  
% nStart (optional) is the nStart passed to
%   SmithCapacity.computeCapacityOnlyAmplitudeConFast. (Default: 2)
%
% The outputs are computed element-wise for each value in A.

if(~issorted(A))
    error('Parameter A must be sorted so this method can work effiently.');
end

% preallocate output
C = zeros(size(A));
n = zeros(size(A));
t = zeros(size(A));

if(nargin < 2)
    nStart = 2;
end

for ii = 1:length(A)
    
    tic;
    [C(ii), poi, ~] = ...
        SmithCapacity.computeCapacityOnlyAmplitudeConFast(A(ii), nStart);
    t(ii) = toc;
    
    % Update nStart.
    % If this value of A(ii) has n=n(ii) as the optimal n, then later
    % (larger) values of A(jj), where A(jj)>A(ii), will have an optimal
    % n>=n(ii).
    n(ii) = length(poi);
    if(n(ii) > nStart)
        nStart = n(ii);
    end
    
end % for ii

end


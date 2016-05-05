function [C, n] = computeTableOfCapacity(A)
% computeTableOfCapacity computes a table of capacity as a function of A.
%
%   [C, n] = computeTableOfCapacity(A)
%
% C is the capacity as computed via Smith1971.
% n is the optimum number of discrete constellation points in X.
%
% A (vector) is the amplitude limit of X (in Y=X+N).  A should be sorted in
%   ascending order.  
%
% The outputs are computed element-wise for each value in A.

if(~issorted(A))
    error('Parameter A must be sorted so this method can work effiently.');
end

% preallocate output
C = zeros(size(A));
n = zeros(size(A));

nStart = 2;

for ii = 1:length(A)
    
    [C(ii), poi, ~] = ...
        SmithCapacity.computeCapacityOnlyAmplitudeConFast(A(ii), nStart);
    
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


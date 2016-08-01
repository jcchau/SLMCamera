function Q = computeQTTransform(G)
% computeQTTransform computes the matrix Q for the Q^T transform using the
% procedure described in p.122-124,127 in lab book #4 (i.e., via
% Gram-Schmidt orthonormalization, but skipping columns of G that are
% linear combinations of previous columns of G).  
%
% Note that this method produces only generates rank(G) columns of Q.  
%
%   Q = computeQTTransform(G)
%
% G is the channel matrix for which Q is generated.  
% Q is the matrix whose transpose is used in the Q^T transform.  
%
% I recommend simplifying G using MIMOCapacity.simplyChannelMatrix before
% passing G to this function.  Although the Q' Transform should also get
% rid of zero rows and columns, simplifyChannelMatrix can combine
% transmitters in a way that this method can't.  
%
% Special case: if G is all zeros, Q will have zero columns.  

[n_r, n_t] = size(G);

% Preallocate matrix Q.
% Will need to remove extra columns later if G does not have full column
% rank.
Q = zeros(n_r, n_t);

% Number of columns computed for Q.
num_col_Q = 0;

for iG = 1:n_t
    colG = G(:, iG);
    
    % The original magnitude of colG before orthogonalization.
    orig_magnitude = norm(colG);
    
    % Orthogonalize colG with respect to the previous columns.  
    for iQ = 1:num_col_Q
        colG = colG - dot(colG, Q(:,iQ)) * Q(:,iQ);
    end % for iQ
    
    % This non-zero threshold is chosen to avoid mistaking colG as linearly
    % independent due to a floating-point math precision error.  
    if(norm(colG) > MIMOCapacity.MinOrthogonalPartFactor * orig_magnitude)
        % Normalize the orthogonalized colG and store it as a new column of
        % Q.  
        num_col_Q = num_col_Q + 1;
        Q(:,num_col_Q) = colG ./ norm(colG);
    end % if norm(colG)
end % for iG

% Resize Q to only keep the num_col_Q columns that were computed.  
Q = Q(:,1:num_col_Q);

end


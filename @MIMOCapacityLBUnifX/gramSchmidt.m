function W = gramSchmidt(V)
% gramSchmidt performs Gram-Schmidt orthogonalization on the columns of V.
%
%   W = gramSchmidt(V)

[numrows, numcols] = size(V);

if(rank(V) < numcols)
    % If V does not have full column rank, the remaining columns would
    % likely not make sense.  
    error('V does not have full column rank.');
end

W = zeros(numrows, numcols);

for iV = 1:numcols
    w_tilde = V(:,iV);
    
    % Orthogonalize w_tilde.
    for iW = 1:iV-1
        w_tilde = w_tilde - dot(w_tilde, W(:,iW)) * W(:,iW);
    end % for iW
    
    % Normalize w_tilde and save as W(:,iV).
    W(:,iV) = w_tilde ./ norm(w_tilde);
end % for iV

end


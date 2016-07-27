function [Q, G_B] = simplifyChannelMatrix(G)
% simplifyChannelMatrix applies the process in lab book #4, p.121-127 to
% reduce the number of dimensions in the channel matrix G.
%
%   [Q, G_B] = MIMOCapacity.simplifyChannelMatrix(G)
%
% Q is the matrix that G_B is multiplied by in the "Q^T transform" to
%   reduce the number of dimensions used to represent the received signals
%   to the rank of the channel matrix G_B.  
% G_B is the simplified channel matrix after step B on p.121 of lab book
%   #4.  It is the channel matrix G with zero rows and columns removed and
%   with the columns that are multiples of each other combined.  
%
% G is the original, unsimplified channel matrix.  
%
% This simplification does not change the capacity of the channel when
% - the noise w is isotropic (i.e., the noise for each receiver element is
%   i.i.d. AWGN) and
% - the transmitter elements are identical.  

%% Step A: remove zero rows and columns
G_A = MIMOCapacity.removeZeroRowsAndCols(G);

%% Step B: combine columns that are multiples of each other

[n_r, n_t] = size(G_A);

% Preallocate G_B.
% We will remove extra columns later.
G_B = zeros(n_r, n_t);
num_cols_GB = 0;

% Preallocate matrix for normalized G_B.
% To avoid needing to repeatedly normalize the columns of G_B.  
normalized_GB = zeros(n_r, n_t);

for ii = 1:n_t

    % The minimum gain for a component of this column that would be
    % considered non-negligible.  
    min_component_gain = norm(G_A(:,ii)) / 1e4;
    
    combined = false;
    
    for jj = 1:num_cols_GB
        % Determine the part of column ii of G_A (as a vector) that is not
        % parallel to column jj of G_B.  
        indep_component = G_A(:,ii) - ...
            dot(G_A(:,ii), normalized_GB(:,jj)) .* normalized_GB(:,jj);
        
        if(norm(indep_component) < min_component_gain)
            % Only a negligible component of the column ii of G_A is
            % independent of column jj of G_B.  
            
            % Instead of making this column in G_A a separate column in
            % G_B, add this column of G_A to the column of G_B that this
            % column is a multiple of.  
            G_B(:,jj) = G_B(:,jj) + G_A(:,ii);
            normalized_GB(:,jj) = G_B(:,jj) ./ norm(G_B(:,jj));
            
            combined = true;
            break;
        end % if(norm(indep_component) < min_component_gain)
    end % for jj = 1:num_cols_GB
    
    if(~combined)
        num_cols_GB = num_cols_GB + 1;
        
        G_B(:, num_cols_GB) = G_A(:, ii);
        normalized_GB(:, num_cols_GB) = G_B(:, num_cols_GB) ./ ...
            norm(G_B(:, num_cols_GB));
    end % if(~combined)
    
end % for ii = 1:n_t

% Keep only the columns of G_B that we filled in
G_B = G_B(:, 1:num_cols_GB);

clear n_r n_t normalized_GB

%% Step C: compute Q for the Q^T transform

Q = MIMOCapacity.computeQTTransform(G_B);

% TODO: Test this method and computeQTTransform.

end


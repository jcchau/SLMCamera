function [umin, umax] = computeUExtremes(F, xmax)
% computeUExtremes computes the extreme values of u=F*x along each
% dimension of u where each element of column vector x is in the range
% [0,1] and where F may contain negative values.  
%
% This method is useful to determine the bin boundaries for the received
% signal z in the channel z=Q'*y=Q'*(G*x+w), where u would be the
% noise-free signal component of z.
% Although every element of channel matrix G is positive in an IM/DD
% optical channel, the channel matrix after the Q'-transform may have
% negative elements.
% As a result, components of u may be negative. 
%
%   [umin, umax] = MIMOCapacity.computeUExtremes(F)
%
% F is the channel matrix after the Q'-transform.  F=Q'*G, where G is the
%   channel matrix before the Q'-transform.  
% xmax is the maximum value of x.  If xmax is a vector (i.e., if the
%   maximum transmit signal needs to be individually specified for each
%   transmitter, then xmax must have n_t elements.  xmax may be either a
%   scalar, a row vector, or a column vector.  
%
% umin is a n_r-element column vector of the minimum value of u.
% umax is a n_r-element column vector of the maximum value of u.  
%
% This algorithm is described in lab book #4, p. 129-133.

[n_r, ~] = size(F,1);

% Take into account xmax. 
% This step allows us to normalize the range of x for the rest of this
% method: 0<=x(i)<=1 for each element of x.  
F = F * diag(xmax);

% Calculate umin and umax
umin = zeros(n_r,1);
umax = zeros(n_r,1);
for d = 1:n_r
    
    % To calculate umin for dimension d, add together each of the negative
    % elements in row d of F.  
    umin(d) = sum(F(d, F(d,:)<0));
    
    % And to calculate umax for dimension d, add together each of the
    % positive elements in row d of F.  
    umax(d) = sum(F(d, F(d,:)>0));
end % for d

end


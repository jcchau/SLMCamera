function [C, poi, voi] = computeCapacityOnlyAmplitudeCon(Alim, delta)
% computeCapacityOnlyAmplitudeCon computes the capacity of scalar channel
% Y = X + N, where noise N is assumed to be Gaussian with zero mean and
% unit variance, and X is constrained to values in [-A, A].
%
% This method uses the algorithm presented in Smith1971.  No constraint on
% variance is assumed.
%
%   [C, poi, voi] = computeCapacityOnlyAmplitudeCon(Alim, delta)
%
% C is the calculated capacity in nats.
% poi (vector) is the points of increase in F_o, the optimal distribution
%   of X. 
% voi (vector) is the probability of each point of increase in poi.
%
% Alim is the amplitude limit of X, so X is constrained to [-A, A].
% delta is the small value by which the algorithm increments A until A
%   reaches Alim.  
%
%% Channel normalization:
% From Smith1969 p. 11--12.
%
% A is defined as (b-a)/(2*sigma_N) for X in [a,b], and N with mean mu and
% variance sigma_N^2.
%
% X' = (X - (a+b)/2) / sigma_N
% N' = (N - mu) / sigma_N
% Y' = X' + N'
%
% Then
% H(Y') = H(Y) - log(sigma_N^2)
% H(N') = H(N) - log(sigma_N^2)
% and
% I(X;Y) = I(X';Y')

% Start with A<=1.6, because then, from Smith1969 (p.51), we know that the
% optimium probability distribution function (or CDF) is F_o(x) = 0.5 *
% U(X+A) + 0.5 * U(X-A) where U is the unit step function.  
A = min(Alim, 1.6);
n = 2;

% Define the distribution F with poi (calculated in the while loop) and
% voi.
% Value of increase (in F_o) at each point of increase (for the two points)
voi = [0.5, 0.5]';

% For n==2.
% The optimal distribution F_o for n=2 is much easier to compute; the code
% in this while loop is optimized for the n=2 case.  
while(true)
    % point of increase
    poi = [-A, A]';
    
    % Mutual information for the optimal distribution of X (given n=2).
    % poi and voi define the optimal distribution for n==2 (farthest
    % equally-probably points of increase).
    I_Fo = SmithCapacity.I(poi, voi);
    
    if(SmithCapacity.checkCorollary1(A, poi, voi, I_Fo))
        if(A == Alim)
            % Found the optimal (maximum) mutual information for A=Alim.
            % Return this mutual information as the capacity.
            C = I_Fo;
            return;
        else
            % Increment A by delta
            A = min(A+delta, Alim);
        end
    else
        % n==2 is no longer optimal, increment n and move on to the next
        % loop (which handles cases where n>2).
        n = n+1;
        break;
    end
end % while (the n=2 loop)

% For n>2
while(true)
    
    %% Compute the optimal Fo given n and A.  
    % Maximize I(Z) by minimizing -I(Z) using fmincon.
    
    % Constraints on Z (lab book #4, p.85-86.)
    Aeq = [ones(1,n), zeros(1,n)];
    beq = 1;
    cA1 = -eye(n, 2*n);
    cA2 = [zeros(n), eye(n)];
    cA3 = [zeros(n), -eye(n)];
    cA = [cA1; cA2; cA3];
    b = [zeros(n,1); repmat(A,2*n,1)];
    
    % A starting point for the optimization: equally-spaced and
    % equally-likely points of increase, spread out from -A to A.
    Z_init = [repmat(1/n, n, 1); (linspace(-A, A, n))'];
    
    % The optimal Z that minimizes -I(Z).
    [Zo, negI_Fo, exitflag] = ...
        fmincon(@(z) -SmithCapacity.I_Z(z), Z_init, cA, b, Aeq, beq);
    if(exitflag<1)
        error('Optimization to find F_o failed.');
    end
    I_Fo = -negI_Fo;
    voi = Zo(1:n);
    poi = Zo(n+1:end);
    
    %% Check for optimality (Smith1971 Corollary 1)
    % Note that by default, for the fmincon function used above:
    % Step Tolerance: 1e-10;
    % Function Tolerance: 1e-6;
    % Optimality Tolerance: 1e-6.  
    % SmithCapacity.checkCorollary1 should allow for these tolerances.
    if(SmithCapacity.checkCorollary1(A, poi, voi, I_Fo))
        if(A == Alim)
            C = I_Fo;
            return;
        else
            A = min(A+delta, Alim);
        end
    else
        n = n+1;
    end
end % while (the n>2 loop)

end

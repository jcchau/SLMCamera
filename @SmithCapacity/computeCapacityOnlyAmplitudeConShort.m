function [C, poi, voi] = computeCapacityOnlyAmplitudeConShort(Alim)
% computeCapacityOnlyAmplitudeConShort computes the capacity of scalar
% channel Y = X + N, where noise N is assumed to be Gaussian with zero mean
% and unit variance, and X is constrained to values in [-A, A].
%
% This method should yield the same result as
% computeCapacityOnlyAmplitudeCon, but skips straight to A=Alim, which
% should save a lot of time for large values of Alim.  
%
% This method uses the algorithm presented in Smith1971.  No constraint on
% variance is assumed.
%
%   [C, poi, voi] = computeCapacityOnlyAmplitudeCon(Alim)
%
% C is the calculated capacity in nats.
% poi (vector) is the points of increase in F_o, the optimal distribution
%   of X. 
% voi (vector) is the probability of each point of increase in poi.
%
% Alim is the amplitude limit of X, so X is constrained to [-A, A].
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

A = Alim;
n = 2;

%% For n==2.
% The optimal distribution F_o for n=2 is much easier to compute; the code
% in this section is optimized for the n=2 case.  

% Define the distribution F with point of increase (poi) and value of
% increase (voi).
poi = [-A, A]';
voi = [0.5, 0.5]';

% Mutual information for the optimal distribution of X (given n=2).
% poi and voi define the optimal distribution for n==2 (farthest
% equally-probably points of increase).
I_Fo = SmithCapacity.I(poi, voi);

if(SmithCapacity.checkCorollary1(A, poi, voi, I_Fo))
    C = I_Fo;
    return;
end

% n==2 is no longer optimal, increment n and move on to the next
% loop (which handles cases where n>2).
n = n+1;

%% Modify the optimization options.  
% The termination of the optimization is determined by TolX (the
% termination tolerance on x) and TolFun (the termination tolerance on the
% function).  
%
% 1) From my past experience, I think we can get Z in I(Z) to a precision
% of 1e-12 (without running into precision problems with MATLAB's
% double-precision floating point math).  So set TolX to 1e-12.
% 2) For this much precision in Z to be meaningful, also set TolCon to
% 1e-12 (otherwise, getting the optimal Z to within 1e-12 wouldn't mean
% much if the total probability in voi exceeds 1 by more than 1e-12).
% 3) checkCorollary1 expects I(F_o) and i(x;F) to be calculated to each be
% calculated to (at worst) \pm 0.5e-6 precision.  However, I don't really
% want to compromise the precision of 1) and 2) by terminating the
% optimization when get this precision on I_Fo for a particular A and n.
% So, also set TolFun to 1e-12 so that TolX is what primarily determines
% when the optimization ends.
% 
% Decrease TolX from 1e-10 to 1e-15. 
% Decrease TolCon from 1e-6 to 1e-12
% Decrease TolFun from 1e-6 to 1e-12.
%
% Note though, that even with TolX = 1e-12 being the terminating condition
% for the optimization, the optimal poi found was still found to be 1.3e-7
% off from the true optimal poi.  (Determined by checking poi(1) and
% poi(end) against -A and A.)  And the middle points of increase are also
% off (determined by checking symmetry); decreasing TolX to 1e-15 seems to
% improve this precision slightly.  Decreasing TolX further to 1e-16 seems
% to worsen precision (probably because fmincon can't recognize when it
% reaches the optimal and due to precision errors in MATLAB's
% floating-point math).
ooptions = optimoptions('fmincon', ...
    'TolX', 1e-15, ...
    'TolFun', 1e-12, ...
    'TolCon', 1e-12, ...
    'MaxFunEvals', 1e5, ... % Need to increase MaxFunEvals
    'MaxIter', 1e4, ...
    'UseParallel', true, ...
    'Display', 'off'); 

%% For n>2
while(true)
    %% Compute the optimal Fo given n and A.  
    % Maximize I(Z) by minimizing -I(Z) using fmincon.

    fprintf('n=%d ', n);
    
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
    
    fprintf('fmincon... ');
    
    % The optimal Z that minimizes -I(Z).
    [Zo, ~, exitflag] = ...
        fmincon(@(z) -SmithCapacity.I_Z(z), Z_init, cA, b, Aeq, beq, ...
        [], [], [], ooptions);
    if(exitflag<1)
        error('Optimization to find F_o failed. A=%.1f, n=%d.', ...
            A, n);
    end
    voi = Zo(1:n);
    poi = Zo(n+1:end);
    
    % Fix precision error for the farthest points of increase in the
    % optimization result.
    poi(1) = -A;
    poi(end) = A;
    I_Fo = SmithCapacity.I(poi, voi);
    
    fprintf('checking...\n');
    
    %% Check for optimality (Smith1971 Corollary 1)
    % Note that SmithCapacity.checkCorollary1 should allow for tolerances.
    if(SmithCapacity.checkCorollary1(A, poi, voi, I_Fo))
        C = I_Fo;
        return;
    else
        % The current value of n is not optimal for A.
        
        % Catch runaway n (due to programming or precision error).
        if(n > 2*A)
            error('n=%d is far too large for A=%.1f.', n, A);
        elseif(n > 1.4*A+1)
            warning('n=%d seems to be too large for A=%.1f.', n, A);
        end
        
        % increment n
        n = n+1;
        
        % Update ooptions.TypicalX to set feature scaling in fmincon.
        ooptions.TypicalX = [repmat(1/n, n, 1); repmat(n/2, n, 1)];
    end % if-else (SmithCapacity.checkCorollary1)
    
end % while (the n>2 loop)

end

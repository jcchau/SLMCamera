function optimal = checkCorollary1(A, poi, voi, I_Fo)
% ??? How do we check that i(x, poi, voi) <= I_Fo for all x in
% [-A, A]?
%
% poi is assumed to be ascending order

% Given that i(x;F) = h(x;F) - D.
% Local maxima of i(x,F) (with respect to x) are the same as
% the local maxima of h(x;F).
% h(x;F) = - integral(p_N(y-x)*log(p_Y(y;F)), -Inf, Inf).
%
% This is equal to the convolution of p_N(y) against
% -log(p_Y(y;F)).  
%
% Thus, all local maxima of h(x;F) are also local maxima of
% -log(p_Y(y;F)).  
% Every local maxima of -log(p_Y(y;F)) is a local minima of
% log(p_Y(y;F)), and hence, a local minima of p_Y(y;F).  
%
% p_Y(y;F) is a weighted sum of Gaussians with mean poi.  
% If two adjacent points of increase are far apart, then in
% between those two points of increase, a single local minimum
% exists.  
%
% If two adjacent points of increase are very close, then the
% peaks of the Gaussian PDF for each point of increase merge to
% form a new peak in between the points of increase.  In this
% case, the local minima between the two points of increase are
% at the points of increase.  
%
% Finally, if two adjacent points of increase are moderately
% close, two peaks exist between the points of increase.  Thus,
% there are three local minima: one at each point of increase,
% and one in between the two peaks.  
%
% Use MATLAB's fminbnd to find the local minima of
% log(p_y(y;F)). 

% Due to precision errors, MATLAB may calculate i(x;F)>I(F_o) even if
% i(x;F) is actually equal or slightly less than I(F_o).  To avoid this
% error, require i(x;F) to exceed I(F_o) by at least TOL_I before deciding
% that i(x;F)>I(F_o).  
TOL_I = 1e-9;
% As a shortcut, just add TOL_I to I_Fo in this function.
I_Fo = I_Fo + TOL_I;

% In addition, the optimization method fmincon will not return precisely
% the correct poi and voi.  
% Among the default tolerances, are:
% - Step Tolerance (where the step is of Z for function I(Z)): 1e-10,
% - Function Tolerance (where the function is I(F)): 1e-6; and
% - Optimality Tolerance: 1e-6.  
% Not understanding how optimality tolerance affects this check, ignore it
% for now.  
% As a loose ballpark estimate of the tolerance, allow i(x;F) to exceed
% I(F_o) by two times the function tolerance.  
I_Fo = I_Fo + 2e-6;
% And allow a tolerance of TOL_Ext when checking that the extreme values of
% poi are -A and A.  
TOL_Ext = 2e-6;

optimal = true;

% First check the extremes.
% 
if(abs(-A-poi(1))>TOL_Ext || abs(A-poi(end))>TOL_Ext)
    % Know that in the optimal distribution, we have a POI at
    % each of the extremes.
    optimal = false;
elseif(any(SmithCapacity.i(poi, poi, voi) > I_Fo))
    % For if the points of increase are very close or
    % moderately close, check at these points of increase if
    % Corollary 1 from Smith 1971 is violated.  
    % Don't need to separately check i(x;F) at x=A and x=-A
    % because we've verified that they're at a point of
    % increase.  
    optimal = false;
else 
    for ii = 1:length(poi)-1
        x1 = poi(ii);
        x2 = poi(ii+1);

        % The local minimum if poi are sufficiently far apart.
        localmaxat = fminbnd( ...
            @(x) SmithCapacity.p_Y(x, poi, voi), ...
            x1, x2);

        if(SmithCapacity.i(localmaxat, poi, voi) > I_Fo)
            optimal = false;
            break
        end

        %% Moderately close case
        % If localmaxat equals x1 or x2, then we may have the
        % case where these two points of increase are
        % moderately close, where there may be another local
        % minimum that we haven't checked yet, hiding between
        % two local maxima.  

        % From fminbnd's documentation:
        % If the minimum actually occurs at x1 or x2, fminbnd
        % returns a point x in the interior of the interval
        % (x1, x2) that is close to the minimizer. In this
        % case, the distance of x from the minimizer is no more
        % than 2*(TolX + 3*abs(x)*sqrt(eps)).
        % The default TolX is 1e-4.
        dist_threshold = ...
            2*(1e-4 + 3*abs(localmaxat) * sqrt(eps));
        if(localmaxat-x1 < dist_threshold || ...
                x2-localmaxat < dist_threshold)
            % I don't know a good way to find the hidden
            % minimum between the maxima, so just sample 1000
            % times between x1 and x2.
            % TODO: find a better way.
            % Perhaps gradenient ascent up to the peaks from
            % both ends, and then find the minimum.  

            x = linspace(x1, x2, 1000);
            iSC = SmithCapacity.i(x, poi, voi);
            if(any(iSC > I_Fo))
                optimal = false;
            end

        end % if is within dist_threshold of x1 or x2
    end % for ii
end % if-else

end % function checkCorollary1

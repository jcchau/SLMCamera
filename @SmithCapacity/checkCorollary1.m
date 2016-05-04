function optimal = checkCorollary1(A, poi, voi, I_Fo)
% ??? How do we check that i(x, poi, voi) <= I_Fo for all x in
% [-A, A]?
%
% poi is assumed to be ascending order

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
% - Step Tolerance (where the step is of Z for function I(Z)),
% - Function Tolerance (where the function is I(F)), and
% - Optimality Tolerance (in newer versions of MATLAB).  
% Not understanding how optimality tolerance affects this check, ignore it
% for now.  
% As a loose ballpark estimate of the tolerance, allow i(x;F) to exceed
% I(F_o) by a little.
I_Fo = I_Fo + 1e-6;

% For any A>0, optimal n is at least 2.
if(length(poi)<2 && A>0)
    optimal = false;
    return
end

% Given that i(x;F) is fairly smooth with respect to x, I think we should
% probably catch any i(x;F) > I_Fo if we just sample between -A and A 100
% times per point of increase.
% TODO: There's probably a better way to do this (e.g., finding all of the
% maxima of i(x;F) and just checking there), but I haven't figured out how
% to reliably do that yet.  

nsamples = 100*numel(poi);

optimal = all( ...
    SmithCapacity.i(linspace(-A, A, nsamples), poi, voi) <= I_Fo);

end % function checkCorollary1

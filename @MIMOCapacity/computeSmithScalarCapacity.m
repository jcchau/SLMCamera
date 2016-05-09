function C = computeSmithScalarCapacity(x_max, sigma_w)
% computeSmithScalarCapacity computes the capacity of optical scalar
% channel y=x+w, where x is in [0,x_max] and w is Gaussian-distributed with
% zero mean and variance sigma_w^2.  
%
%   C = computeSmithScalarCapacity(x_max, sigma_w)
%
% C is the computed capacity in nats.
%
% x_max is the upper bound of possible values of x.
% sigma_w is the standard deviation of the Gaussian noise w.
%
% x_max and sigma_w may be vectors of the same size; in this case, C will
% be computed element-wise for each (x_max, sigma_w)-pair.  

Alim = x_max ./ (2 * sigma_w);

try
% Use the fast algorithm (instead of the one that matches the flow-chart
% from Smith1971).  
C = SmithCapacity.computeCapacityOnlyAmplitudeConFast(Alim);
catch me
    if(strcmp(me.identifier, 'SmithCapacity:abortA'))
        C = NaN;
    else
        rethrow(me);
    end
end % try-catch

end


function testComputeUExtremesEmpirically(tc)
% Tests method MIMOCapacity.computeUExtremes using a random G and a random
% vector xmax.  
% Evaluates that the extremes, umin and umax hold by trying many extreme
% values of x.

n_r = randi(10);
n_t = randi(10);
G = rand(n_r, n_t);
xmax = rand(n_t,1) ./ rand();

% With a randomly-generated G, MIMOCapacity.simplifyChannelMatrix won't do
% anything, so skip it.  
% Furthermore, MIMOCapacity.simplifyChannelMatrix assumes that xmax is the
% same for every transmitter; we don't use that assumption here.  

Q = MIMOCapacity.computeQTTransform(G);

F = Q' * G;

%% Method under test

[umin, umax] = MIMOCapacity.computeUExtremes(F, xmax);

%% Check that u = F*x is always within [umin(d), umax(d)] for each dim.
% Check for 1000 extreme values of x.

for trial = 1:1000
    % Generate an extreme value of x, where each element is either 0 or the
    % corresponding xmax.  
    x = xmax .* (rand(n_t,1) > 0.5);
    
    u = F * x;
    
    uminOK = (u>=umin) & (u<=umax);
    
    if(~all(uminOK))
        tc.verifyLessThanOrEqual(umin, u, ...
            'umin should be the minimum possible u.');
        tc.verifyGreaterThanOrEqual(umax, u, ...
            'umax should be the maximum possible u.');
        
        % Once we find a failure, no need to continue the test.
        break
    end % if
end % for trial

end


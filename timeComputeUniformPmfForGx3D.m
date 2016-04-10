function timeComputeUniformPmfForGx3D(s_max_nbins, s_threads, s_rseed)
% timeComputeUniformPmfForGx3D times how long it takes to run
% computeUniformPmfForGx on a nbin 3D PMF for a random G.

max_nbins = sscanf(s_max_nbins, '%e');
nthreads = sscanf(s_threads, '%u');

if(nargin >= 3)
    % seed the random number generator if a seed is provided
    rseed = sscanf(s_rseed, '%i');
    rng(rseed);
end

if(~isempty(gcp('nocreate')))
    delete(gcp)
end
parpool(nthreads);
maxNumCompThreads(nthreads);

n_t = 3;
n_r = 3;

nbins = MIMOCapacity.fillMaxNBins(max_nbins, n_r);
fprintf('Total number of bins: %u.\n', prod(nbins));

G = rand(n_r, n_t);
x_max = ones(n_t, 1);
y_min = zeros(1, n_r);
y_max = (G * x_max)';
delta = (y_max - y_min) ./ nbins;

tic

[pmf, reachable] = MIMOCapacity.computeUniformPmfForGx(G, x_max, ...
    y_min, delta, nbins);

toc

if(~sum(pmf(:)) == 1)
    fprintf('pmf does not sum to 1.\n');
end
if(~all(pmf(~reachable) == 0))
    fprintf('Some non-reachable bins have pmf ~= 0.\n');
end

end


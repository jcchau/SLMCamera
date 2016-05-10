function startSimulations(runTrial_cmd, savedir)
% runTrial_cmd is the command used to start the next stage of the
% simulation.  

srand = 20160509;
rng(srand);

ntrials = 100;
srand_trial = randi(intmax, ntrials, 1);

for ii = 1:ntrials
    matlab_cmd = sprintf('%s %d %d "%s"', ...
        runTrial_cmd, ii, srand_trial(ii), savedir);
    job_cmd = sprintf('qsub -N simt%d -cwd -j y -m a -V -b y %s', ...
        ii, matlab_cmd);
    system(job_cmd);
end % for ii

end


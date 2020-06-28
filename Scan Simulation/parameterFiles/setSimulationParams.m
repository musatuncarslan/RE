function [Simparams] = setSimulationParams(MPIparams, Physicsparams)
    Simparams = struct;
    
    fs = Physicsparams.fs;
    ffp_type = MPIparams.ffp_type;
    time = MPIparams.time; % time to traverse whole FOV
   
    
    % find which periods are to be simulated
    if strcmpi(ffp_type, 'linear_rastered')
        simPeriodsz = findPeriodsZ(MPIparams);
        Simparams.simPeriods = simPeriodsz; % periods to simulate
    elseif strcmpi(ffp_type, 'fixed')
    elseif strcmpi(ffp_type, 'complex_rastered')
        simPeriodsx = findPeriodsX(MPIparams);
        simPeriodsz = findPeriodsZ(MPIparams);
        Simparams.simPeriods = intersect(simPeriodsx,simPeriodsz,'stable'); % periods to simulate
    end
    
    % adjust the simulation sampling frequency so that each period has
    % integer number of samples and Physics (also simulation) sampling 
    % frequency is divisible by both MPI sampling frequency and drive
    % field frequency.
    f_drive = MPIparams.f_drive;
    fs_mpi = MPIparams.fs;
    fs_mpi = fs_mpi-mod(fs_mpi, f_drive);
    fs_phsy = Physicsparams.fs;
    fs_phsy = fs_phsy + fs_mpi/2;
    fs_phsy = fs_phsy - mod(fs_phsy, fs_mpi);
    Simparams.fs_mpi = fs_mpi;
    Simparams.fs_phsy = fs_phsy;
    Simparams.downsample = fs_phsy/fs_mpi;
    
    numPeriods = time*f_drive; % number of drive periods during whole FOV scan (this is hopefully an integer)
    if (floor(numPeriods)~=numPeriods)
        error('Error. Make sure that time*f_drive is an integer. \n\n \t You can slightly adjust "%s" for this.', 'time')
    end
    Simparams.samplePerPeriod = fs*time/numPeriods; % number of samples per period
    
    
    
%     
%     ffp_type = MPIparams.ffp_type;
%     
%     if strcmpi(ffp_type, 'linear_rastered')
%         [startIter, endIter, samplesPerIter] = linearRasteredIteration(MPIparams, Physicsparams);
%         numIters = endIter - startIter+1;
%         Simparams.startIter = startIter;
%         Simparams.endIter = endIter;
%     elseif strcmpi(ffp_type, 'fixed')
%         numIters = MPIparams.cycle;
%         total_time = numIters/MPIparams.f_drive;
%         samplesPerIter = Physicsparams.fs*total_time/numIters;
%         Simparams.startIter = 1;
%         Simparams.endIter = Simparams.startIter + numIters - 1;
%     elseif strcmpi(ffp_type, 'triangular')
%     end
%     
% 
%     Simparams.numSamplesPerIter = samplesPerIter;
%     Simparams.numIters = numIters;
% 
%     Simparams.downsample = Physicsparams.fs/MPIparams.fs; % downsample ratio

end
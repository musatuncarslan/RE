function [Simparams, MPIparams] = setSimulationParams(MPIparams, Physicsparams)
    Simparams = struct;
    
    ffp_type = MPIparams.ffp_type;
    
    
    % find which periods are to be simulated
    if strcmpi(ffp_type, 'linear_rastered')
        MPIparams = adjustFOVz(MPIparams); % adjust z-FOV if necessary
        zPeriodsSim = simZperiods(MPIparams);
        Simparams.simPeriods = zPeriodsSim; % periods to simulate
    elseif strcmpi(ffp_type, 'fixed')
    elseif strcmpi(ffp_type, 'complex_rastered')
        MPIparams = adjustFOVz(MPIparams); % adjust z-FOV if necessary
        zPeriodsSim = simZperiods(MPIparams);
        xPeriodsSim = simXperiods(MPIparams);
        Simparams.simPeriods = intersect(xPeriodsSim,zPeriodsSim,'stable'); % periods to simulate
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
    
    MPIparams.fs = fs_mpi;
    
    time = MPIparams.time; fs = Simparams.fs_phsy;
    Simparams.samplePerPeriod = fs_phsy/f_drive; % number of samples per period
    
    [Simparams.div, Simparams.divNum] = partitionPeriods(Simparams);



end
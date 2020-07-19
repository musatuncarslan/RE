function [Simparams, MPIparams] = setSimulationParams(MPIparams, Physicsparams)
    Simparams = struct;
    
    
    if ((MPIparams.Rs(3) ~= 0) && (MPIparams.Rs(1) == 0))
        MPIparams = adjustFOVz(MPIparams); % adjust z-FOV if necessary
        zPeriodsSim = simZperiods(MPIparams);
        Simparams.simPeriods = zPeriodsSim; % periods to simulate
    elseif ((MPIparams.Rs(3) ~= 0) && (MPIparams.Rs(1) ~= 0))
        MPIparams = adjustFOV(MPIparams); % adjust z-FOV if necessary
        
        MPIparams.traversedFOVz = [0.00 MPIparams.FOV_z/(MPIparams.numTrianglePeriods*2)];
        MPIparams.traversedFOVx = [-0.01 0.01];
        
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
    
    Simparams.samplePerPeriod = fs_phsy/f_drive; % number of samples per period
    
    [Simparams.div, Simparams.divNum] = partitionPeriods(Simparams);



end
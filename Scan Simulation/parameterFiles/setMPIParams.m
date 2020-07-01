function MPIparams = setMPIParams(Physicsparams, ffp_type, rs)
    MPIparams = struct;

    MPIparams.ffp_type = ffp_type;
    if strcmp(MPIparams.ffp_type, 'fixed') == 1
        MPIparams.slewRate = 0;
    elseif strcmp(MPIparams.ffp_type, 'linear_rastered') == 1
        MPIparams.slewRate = rs; % slew rate (T/s)
    elseif strcmp(MPIparams.ffp_type, 'complex_rastered') == 1
        MPIparams.slewRate = rs; % slew rate (T/s)
    end
    % gradients (T/m) (current scanner)
    MPIparams.Gxx = 4.8;
    MPIparams.Gyy = 2.4;
    MPIparams.Gzz = 2.4;
    
    MPIparams.Bp = 10e-3; % Drive field (T)
    MPIparams.f_drive = 20e3; % drive field frequency

    
    Hp=MPIparams.Bp/Physicsparams.mu0; % magnetization moment
    G=MPIparams.Gzz/Physicsparams.mu0; % gradient
    MPIparams.driveMag=Hp/G; % extent of the drive field

    MPIparams.fs = 2e6; % sample frequency of the MPI system (Hz)
    MPIparams.FOV_z = 0.06; % FOV in z-axis (meters) (bore axis)
    MPIparams.FOV_x = 0.05; % FOV in x-axis (meters)
    
    if strcmp(MPIparams.ffp_type, 'linear_rastered') == 1
        % for linear rastered
        MPIparams.time = MPIparams.FOV_z*MPIparams.Gzz*(1/MPIparams.slewRate); % time (seconds)
        MPIparams.traversedFOVz = [-0.01 0.01]; % traversed fov in the simulation in z-axis (m)
    elseif strcmp(MPIparams.ffp_type, 'complex_rastered') == 1
        % for complex rastered
        MPIparams.time = MPIparams.FOV_z*MPIparams.Gzz*(1/MPIparams.slewRate); % time (seconds)
        MPIparams.numTrianglePeriods = 5;
        MPIparams.traversedFOVz = [-0.00 0.006];
        MPIparams.traversedFOVx = [-0.02 0.02];
    end

    

    
    % related to time constant estimation
    MPIparams.interp_coeff = 10;
   
end
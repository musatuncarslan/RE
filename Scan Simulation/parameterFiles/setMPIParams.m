function MPIparams = setMPIParams(Physicsparams, rs)

    MPIparams.Rs = rs; % slew rate vector [x,y,z] (T/s)
    
    % gradients (T/m) (current scanner)
    MPIparams.Gxx = 4.8;
    MPIparams.Gyy = 2.4;
    MPIparams.Gzz = 2.4;
    
    MPIparams.Bp = 10e-3; % Drive field (T)
    MPIparams.f_drive = 10e3; % drive field frequency

    
    Hp=MPIparams.Bp/Physicsparams.mu0; % magnetization moment
    G=MPIparams.Gzz/Physicsparams.mu0; % gradient
    MPIparams.driveMag=Hp/G; % extent of the drive field

    MPIparams.fs = 2e6; % sample frequency of the MPI system (Hz)
    MPIparams.FOV_z = 0.06; % FOV in z-axis (meters) (bore axis)
    MPIparams.FOV_x = 0.05; % FOV in x-axis (meters)
   
    MPIparams.time = MPIparams.FOV_z*MPIparams.Gzz*(1/MPIparams.Rs(3)); % time (seconds)
    MPIparams.numTrianglePeriods = 5;


end


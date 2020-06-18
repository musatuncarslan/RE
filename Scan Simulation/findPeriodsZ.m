function simPeriodsz = findPeriodsZ(MPIparams)

    % calculation of necessary periods in the z-direction for partial FOV
    % simulation
    f_drive = MPIparams.f_drive;
    time = MPIparams.time; % time to traverse whole FOV
    numPeriods = time*f_drive; % number of drive periods during whole FOV scan (this is hopefully an integer)
    if (floor(numPeriods)~=numPeriods)
        error('Error. Make sure that time*f_drive is an integer. \n\n \t You can slightly adjust "%s" for this.', 'time')
    end
    FOV_z = MPIparams.FOV_z;
    robotSpeed_z = FOV_z/time; % robot arm movement speed in z direction (m/s)
    traversedFOVz = [-0.01 0.01]; % z-axis FOV that is of interest (m)

    zi = robotSpeed_z*time/numPeriods;
    z11 = round((traversedFOVz(1)+FOV_z/2)/zi);
    z12 = round((traversedFOVz(2)+FOV_z/2)/zi)-1;
    simPeriodsz = z11:z12;

end
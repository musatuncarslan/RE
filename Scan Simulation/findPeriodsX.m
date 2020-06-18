function simPeriodsx = findPeriodsX(MPIparams)

    f_drive = MPIparams.f_drive;
    time = MPIparams.time; % time to traverse whole FOV
    numPeriods = time*f_drive; % number of drive periods during whole FOV scan (this is hopefully an integer)
    if (floor(numPeriods)~=numPeriods)
        error('Error. Make sure that time*f_drive is an integer. \n\n \t You can slightly adjust "%s" for this.', 'time')
    end


    % calculation of necessary periods in the x-direction for partial FOV
    % simulation
    time = MPIparams.time; % time to traverse whole FOV
    FOV_x = MPIparams.FOV_x;
    numTrianglePeriods = MPIparams.numTrianglePeriods; % number of triangle zig-zags in whole movement
    numPeriodsPerTriangle = numPeriods/numTrianglePeriods; % number of drive periods per triangle periods 
                                                           % (regarding to the movement in x-direction, 
                                                           % this is also hopefully an even integer, if necessary adjust time so that it is an even integer)
    if (floor(numPeriodsPerTriangle)~=numPeriodsPerTriangle)
        error('Error. Make sure that time*f_drive/numTrianglePeriods is an integer. \n\n \t You can slightly adjust "%s" for this.', 'time')
    end
    if mod(numPeriodsPerTriangle, 2) == 1
        error('Error. Make sure that time*f_drive/numTrianglePeriods is divisible by 2. \n\n \t You can slightly adjust "%s" for this.', 'time')
    end

    robotSpeed_x = FOV_x*numTrianglePeriods*2/time; % robot arm movement speed in x direction (m/s)
    traversedFOVx = [-0.005 0.005]; % x-axis FOV that is of interest (m)

    xi = robotSpeed_x*time/numPeriods;
    x11 = round((traversedFOVx(1)+FOV_x/2)/xi);
    x12 = round((traversedFOVx(2)+FOV_x/2)/xi)-1;
    x21 = round((-traversedFOVx(2)+FOV_x/2)/xi)+numPeriodsPerTriangle/2;
    x22 = round((-traversedFOVx(1)+FOV_x/2)/xi)+numPeriodsPerTriangle/2-1;

    L = ((x12-x11)+(x22-x21))+2;
    simPeriodsx = zeros(1, L*numTrianglePeriods);
    for k=1:numTrianglePeriods
        simPeriodsx((k-1)*L+1:k*L) = ([x11:x12, x21:x22])+(k-1)*numPeriodsPerTriangle; % simulated periods
    end

end



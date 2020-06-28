function [FFPparams] = generateFFP(gpudev, t, MPIparams, Simparams, sigmoidParams)

    FOV_z = MPIparams.FOV_z;
    FOV_x = MPIparams.FOV_x;
    ffp_type = MPIparams.ffp_type;
    time = MPIparams.time;

    driveMag = MPIparams.driveMag;
    f_drive = MPIparams.f_drive;
    

    if strcmp(ffp_type, 'linear_rastered')
        
    elseif strcmpi(ffp_type, 'complex_rastered')
        z = FOV_z/time*t - FOV_z/2; % robot arm movement in z direction w.r.t. time
        p = 1/(MPIparams.numTrianglePeriods/time);
        x = FOV_x*(2*abs(2*(t/p - floor(t/p + 0.5)))-1)/2; % robot arm movement in x direction w.r.t. time
    else
        
    end
    
    drive = cos(2*pi*f_drive*t); % drive field movement
%     drive = sigmf(drive,sigmoidParams)-0.5;
    drive = driveMag*drive/max(drive);
    
    
    FFPparams.drive = drive;
    FFP_z = z+drive; % movement of FFP in z direction (robot arm movement in z direction + drive field movement, since drive field is in z direction as well)
    FFP_x = x;
    
    zDifference = diff(FFP_z, 1, 2);
    xDifference = diff(FFP_x, 1, 2);

    FFP_speed = sqrt(zDifference.^2 + xDifference.^2); 
    rounded = round(atan2d(xDifference, zDifference)*1000)/1000;
    FFP_angle = wrapTo360(rounded);    

    
%     [sigPow, ~] = max(abs(FFP_z));
%     sigma = sigPow/50;
%     FFP_z=FFP_z+sigma*randn(1,length(FFP_z));
    
    FFPparams.FFP_x = FFP_x; wait(gpudev); clear FFP_x;
    FFPparams.FFP_z = FFP_z; wait(gpudev); clear FFP_z;
    FFPparams.FFP_speed = FFP_speed; wait(gpudev); clear FFP_speed;
    FFPparams.FFP_angle = FFP_angle;  
    FFPparams.FFP_uniqueAngle = unique(FFP_angle); wait(gpudev); clear FFP_angle;
    
end
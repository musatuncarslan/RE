function [FFPparams] = generateFFP(gpudev, t, MPIparams, type)

    FOV_z = MPIparams.FOV_z;
    FOV_x = MPIparams.FOV_x;
    time = MPIparams.time;

    driveMag = MPIparams.driveMag;
    f_drive = MPIparams.f_drive;
    
    switch type
        case 'normal'
            if ((MPIparams.Rs(3) ~= 0) && (MPIparams.Rs(1) == 0))
                zSpeed = MPIparams.Rs(3)/MPIparams.Gzz;
                xSpeed = MPIparams.Rs(1)/MPIparams.Gxx;
                z = zSpeed*t - zSpeed*t(end)/2; % robot arm movement in z direction w.r.t. time
                x = xSpeed*t - xSpeed*t(end)/2; % robot arm movement in x direction w.r.t. time
            elseif ((MPIparams.Rs(3) ~= 0) && (MPIparams.Rs(1) ~= 0))
                zSpeed = MPIparams.Rs(3)/MPIparams.Gzz;
                xSpeed = MPIparams.Rs(1)/MPIparams.Gxx;
                z = zSpeed*t - zSpeed*t(end)/2; % robot arm movement in z direction w.r.t. time
                x = xSpeed*t - xSpeed*t(end)/2 + xSpeed/MPIparams.f_drive/2; % robot arm movement in x direction w.r.t. time
            end
        case 'test'
            zSpeed = MPIparams.Rs(3)/MPIparams.Gzz;
            xSpeed = MPIparams.Rs(1)/MPIparams.Gxx;
            z = zSpeed*t - zSpeed*t(end)/2; % robot arm movement in z direction w.r.t. time
            x = xSpeed*t - xSpeed*t(end)/2; % robot arm movement in x direction w.r.t. time
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
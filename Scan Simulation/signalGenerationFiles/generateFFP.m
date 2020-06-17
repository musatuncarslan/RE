function [FFPparams] = generateFFP(t, MPIparams, Simparams, sigmoidParams)

    FOV_z = MPIparams.FOV_z;
    FOV_x = MPIparams.FOV_x;
    ffp_type = MPIparams.ffp_type;
    time = MPIparams.time;

    driveMag = MPIparams.driveMag;
    f_drive = MPIparams.f_drive;
    
    numSamplesPerIter = Simparams.numSamplesPerIter;

    if strcmp(ffp_type, 'linear_rastered')
        FFP_x = repmat(MPIparams.ffpPosition(1), [1, numSamplesPerIter+1+Simparams.downsample*2]); % movement of FFP in x direction, for linear x is constant throughout the line
        z = FOV_z/time*t - FOV_z/2; % robot arm movement in z direction w.r.t. time
    elseif strcmpi(ffp_type, 'fixed')
        numPos = size(MPIparams.ffpPosition, 1);
        FFP_x = zeros(numPos, numSamplesPerIter+1+Simparams.downsample*2);
        z = zeros(numPos, numSamplesPerIter+1+Simparams.downsample*2);
        for k=1:size(MPIparams.ffpPosition, 1)
            FFP_x(k, :) = repmat(MPIparams.ffpPosition(k, 1), [1, numSamplesPerIter+1+Simparams.downsample*2]); % movement of FFP in x direction, for linear x is constant throughout the line
            z(k, :) = repmat(MPIparams.ffpPosition(k, 2), [1, numSamplesPerIter+1+Simparams.downsample*2]);
        end
    else
        
    end
    
    drive = cos(2*pi*f_drive*t); % drive field movement
    drive = sigmf(drive,sigmoidParams)-0.5;
    drive = driveMag*drive/max(drive);
    
    
    FFPparams.drive = drive;
    FFP_z = z+drive; % movement of FFP in z direction (robot arm movement in z direction + drive field movement, since drive field is in z direction as well)
   
    zDifference = diff(FFP_z, 1, 2);
    xDifference = diff(FFP_x, 1, 2);

    FFP_speed = sqrt(zDifference.^2 + xDifference.^2); 
    rounded = round(atan2d(xDifference, zDifference)*100)/100;
    FFP_angle = wrapTo360(rounded);    

    
%     [sigPow, ~] = max(abs(FFP_z));
%     sigma = sigPow/50;
%     FFP_z=FFP_z+sigma*randn(1,length(FFP_z));
    
    FFPparams.FFP_x = FFP_x;
    FFPparams.FFP_z = FFP_z;
    FFPparams.FFP_speed = FFP_speed;
    FFPparams.FFP_angle = FFP_angle;
    FFPparams.FFP_uniqueAngle = unique(FFP_angle);
    
    truncSig = FFPparams.FFP_x(1:Simparams.downsample:end-1);
    FFPparams.FFP_xToFile = gather(truncSig);    
    truncSig = FFPparams.FFP_z(1:Simparams.downsample:end-1);
    FFPparams.FFP_zToFile = gather(truncSig);    
    truncSig = FFPparams.drive(1:Simparams.downsample:end-1);
    FFPparams.FFP_driveToFile = gather(truncSig); 
    truncSig = FFPparams.FFP_speed(1:Simparams.downsample:end);
    FFPparams.FFP_speedToFile = gather(truncSig);
    truncSig = FFPparams.FFP_angle(1:Simparams.downsample:end);
    FFPparams.FFP_angleToFile =  gather(truncSig);
    
    
end
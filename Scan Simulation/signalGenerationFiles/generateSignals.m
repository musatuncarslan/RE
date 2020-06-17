function [signal, SPIOparams] = generateSignals(colinearIMG, transverseIMG, FFPparams, MPIparams, SPIOparams, Simparams)

    signal = struct;

    horizontalPrev = SPIOparams.horizontalPrev;
    verticalPrev = SPIOparams.verticalPrev;
    
    FFP_x = FFPparams.FFP_x;
    FFP_z = FFPparams.FFP_z;
    FFP_speed = FFPparams.FFP_speed;
    FFP_angle = FFPparams.FFP_angle;
    FFP_uniqueAngle = FFPparams.FFP_uniqueAngle;
    
    FOV_x = MPIparams.FOV_x;
    FOV_z = MPIparams.FOV_z;
        
    dx = SPIOparams.dx; 
    dz = SPIOparams.dz;
    r_t = SPIOparams.r_t;

    % generate colinear and transverse signals using interpolation
    x = (-FOV_x/2:dx:FOV_x/2);
    z = (-FOV_z/2:dz:FOV_z/2);
    numAngle = length(FFP_uniqueAngle);
    numSPIO = length(SPIOparams.diameter);
    colinearSignal = gpuArray.zeros(numSPIO, length(FFP_angle));
    transverseSignal = gpuArray.zeros(numSPIO, length(FFP_angle));

    for l=1:numSPIO
        for k=1:numAngle
            angleIdx = (FFP_angle == FFP_uniqueAngle(k));
            colinearSignal(l, angleIdx) = (interp2(x,z,colinearIMG{l, k}, FFP_x(angleIdx), FFP_z(angleIdx), 'linear'));
            transverseSignal(l, angleIdx) = (interp2(x,z,transverseIMG{l, k}, FFP_x(angleIdx), FFP_z(angleIdx), 'linear'));
        end
    end

    % generate horizontal and vertical signals using relaxed colinear and
    % transverse signals.
    horizontalSignal = gpuArray.zeros(size(colinearSignal));
    verticalSignal = gpuArray.zeros(size(colinearSignal));
    for l=1:numSPIO
        for k = 1:numAngle
            angleIdx = (FFP_angle == FFP_uniqueAngle(k));
            horizontalSignal(l, angleIdx) = colinearSignal(l, angleIdx).*FFP_speed(angleIdx).*cosd(FFP_angle(angleIdx)) ...
                - transverseSignal(l, angleIdx).*FFP_speed(angleIdx).*sind(FFP_angle(angleIdx));
            verticalSignal(l, angleIdx) = colinearSignal(l, angleIdx).*FFP_speed(angleIdx).*sind(FFP_angle(angleIdx)) ...
                + transverseSignal(l, angleIdx).*FFP_speed(angleIdx).*cosd(FFP_angle(angleIdx));    
        end
    end

    % create relaxed colinear and transverse signals. Crop the signals 
    % appropriately.
    for k=1:numSPIO
        horizontalSignal_temp = conv([horizontalPrev{k} horizontalSignal(k,:)], r_t{k});
        horizontalPrev{k} = horizontalSignal(k, end-2*Simparams.downsample-length(r_t{k})+1:end-2*Simparams.downsample);
        horizontalSignal(k,:) = horizontalSignal_temp(length(horizontalPrev{k})+1:end-length(r_t{k})+1);
        
        verticalSignal_temp = conv([verticalPrev{k}, verticalSignal(k,:)], r_t{k});
        verticalPrev{k} = verticalSignal(k, end-2*Simparams.downsample-length(r_t{k})+1:end-2*Simparams.downsample);
        verticalSignal(k,:) = verticalSignal_temp(length(verticalPrev{k})+1:end-length(r_t{k})+1);
    end
    

    SPIOparams.horizontalPrev = horizontalPrev;
    SPIOparams.verticalPrev = verticalPrev;
    
    truncSig = colinearSignal(:, 1:Simparams.downsample:end);
    signal.colinearSignal = gather(truncSig);
    truncSig = transverseSignal(:, 1:Simparams.downsample:end);
    signal.transverseSignal = gather(truncSig);
    truncSig = horizontalSignal(:, 1:Simparams.downsample:end);
    signal.horizontalSignal = gather(truncSig);
    truncSig = verticalSignal(:, 1:Simparams.downsample:end);
    signal.verticalSignal = gather(truncSig);
    
end
function [signal, SPIOparams] = generateSe(angleVec, FFPparams, MPIparams, SPIOparams, Simparams, particleNo, div)

    signal = struct;
    
    info = h5info('./temp/IMG.h5');
    dataSize = info.Datasets(1).Dataspace.Size;
    zL = dataSize(1); xL = dataSize(2);
    % divide data into chunks
    numAngle = dataSize(3);
    
    [~,ia] = setdiff(angleVec,FFPparams.FFP_uniqueAngle);
    numAngle = numAngle - length(ia);
    
    numIters = floor(numAngle/div);
    if ((numAngle/div) > numIters)
        numIters = numIters + 1;
    end
    iterVec = [div*ones(1, numIters-1), mod(numAngle, div)];
    idxVec = [0 cumsum(iterVec)];

    horizontalPrev = SPIOparams.horizontalPrev(particleNo);
    verticalPrev = SPIOparams.verticalPrev(particleNo);
    
    FFP_x = FFPparams.FFP_x;
    FFP_z = FFPparams.FFP_z;
    FFP_speed = FFPparams.FFP_speed;
    FFP_angle = FFPparams.FFP_angle;
    FFP_uniqueAngle = FFPparams.FFP_uniqueAngle;
    
    FOV_x = MPIparams.FOV_x;
    FOV_z = MPIparams.FOV_z;
        
    dx = SPIOparams.dx; 
    dz = SPIOparams.dz;
    r_t = SPIOparams.r_t(particleNo);

    % generate colinear and transverse signals using interpolation
    x = gpuArray(-FOV_x/2:dx:FOV_x/2);
    z = gpuArray(-FOV_z/2:dz:FOV_z/2);
    numSPIO = length(SPIOparams.diameter);
    colinearSignal = gpuArray.zeros(numSPIO, length(FFP_angle));
    transverseSignal = gpuArray.zeros(numSPIO, length(FFP_angle));
    horizontalSignal = gpuArray.zeros(size(colinearSignal));
    verticalSignal = gpuArray.zeros(size(colinearSignal));
    
    for l=1:numIters
        colinearIMG = gpuArray(h5read('./temp/IMG.h5','/colinearIMG', [1 1 idxVec(l)+1], [zL xL iterVec(l)]));
        transverseIMG = gpuArray(h5read('./temp/IMG.h5','/transverseIMG', [1 1 idxVec(l)+1], [zL xL iterVec(l)]));
        try
            partialAngle = FFP_uniqueAngle(idxVec(l)+1:idxVec(l+1));
        catch ME
            aa = 5;
        end
        for k=1:length(partialAngle)
            angleIdx = (FFP_angle == partialAngle(k));
            
            % generate relaxed colinear and transverse signals.
            colinearSignal(angleIdx) = (interp2(x,z,colinearIMG(:,:,k), FFP_x(angleIdx), FFP_z(angleIdx), 'linear'));
            transverseSignal(angleIdx) = (interp2(x,z,transverseIMG(:,:,k), FFP_x(angleIdx), FFP_z(angleIdx), 'linear'));

            % generate horizontal and vertical signals using relaxed colinear and
            % transverse signals.
            horizontalSignal(angleIdx) = colinearSignal(angleIdx).*FFP_speed(angleIdx).*cosd(FFP_angle(angleIdx)) ...
                - transverseSignal(angleIdx).*FFP_speed(angleIdx).*sind(FFP_angle(angleIdx));
            verticalSignal(angleIdx) = colinearSignal(angleIdx).*FFP_speed(angleIdx).*sind(FFP_angle(angleIdx)) ...
                + transverseSignal(angleIdx).*FFP_speed(angleIdx).*cosd(FFP_angle(angleIdx));               
        end
    end





    % create relaxed colinear and transverse signals. Crop the signals 
    % appropriately.
    horizontalSignal_temp = conv([horizontalPrev horizontalSignal], r_t);
    horizontalPrev = horizontalSignal(end-2*Simparams.downsample-length(r_t)+1:end-2*Simparams.downsample);
    horizontalSignal = horizontalSignal_temp(length(horizontalPrev)+1:end-length(r_t)+1);

    verticalSignal_temp = conv([verticalPrev, verticalSignal], r_t);
    verticalPrev = verticalSignal(k, end-2*Simparams.downsample-length(r_t)+1:end-2*Simparams.downsample);
    verticalSignal = verticalSignal_temp(length(verticalPrev)+1:end-length(r_t)+1);
    

    SPIOparams.horizontalPrev{particleNo} = horizontalPrev;
    SPIOparams.verticalPrev{particleNo} = verticalPrev;
    
%     signal.colinearSignal = gather(colinearSignal);
%     signal.transverseSignal = gather(transverseSignal);
%     signal.horizontalSignal = gather(horizontalSignal);
%     signal.verticalSignal = gather(verticalSignal);
    
    truncSig = colinearSignal(:, 1:Simparams.downsample:end);
    signal.colinearSignal = gather(truncSig);
    truncSig = transverseSignal(:, 1:Simparams.downsample:end);
    signal.transverseSignal = gather(truncSig);
    truncSig = horizontalSignal(:, 1:Simparams.downsample:end);
    signal.horizontalSignal = gather(truncSig);
    truncSig = verticalSignal(:, 1:Simparams.downsample:end);
    signal.verticalSignal = gather(truncSig);
    
end
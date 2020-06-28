function [signal, SPIOparams] = generateSe(gpudev, FFPparams, MPIparams, SPIOparams, Simparams, particleNo, div)

    signal = struct;

    % initialize parameters
    info = h5info('./temp/PSF.h5');
    dataSize = info.Datasets(2).Dataspace.Size;
    zL = dataSize(1); xL = dataSize(2);
    angleVecSize = info.Datasets(1).Dataspace.Size;
    
    horizontalPrev = SPIOparams.horizontalPrev{particleNo};
    verticalPrev = SPIOparams.verticalPrev{particleNo};
    
    FFP_x = FFPparams.FFP_x;
    FFP_z = FFPparams.FFP_z;
    FFP_speed = FFPparams.FFP_speed;
    FFP_angle = FFPparams.FFP_angle;
    FFP_uniqueAngle = FFPparams.FFP_uniqueAngle;
    
    FOV_x = MPIparams.FOV_x;
    FOV_z = MPIparams.FOV_z;
        
    dx = SPIOparams.dx; 
    dz = SPIOparams.dz;
    r_t = SPIOparams.r_t{particleNo};
    
    
    % divide data into chunks
    angleVec = h5read('./temp/PSF.h5', '/angleVec', [1 1], angleVecSize);
    [~,ia] = intersect(angleVec,gather(FFP_uniqueAngle));
    ia = [1; ia; angleVec(end)]';
    iaP = diff(ia);
    iaPneq1 = iaP(find(iaP ~= 1));
    ianeq1 = ia(find(iaP ~= 1));

    is = ianeq1(1:end-1)+iaPneq1(1:end-1);
    ie = ianeq1(2:end);

    for k=1:length(is)
        tempAngleVec = is(k):ie(k);
        numAngle = length(tempAngleVec);
        numIters{k} = floor(numAngle/div);
        if ((numAngle/div) > numIters{k})
            numIters{k} = numIters{k} + 1;
        end
        iterVec{k} = [div*ones(1, numIters{k}-1), mod(numAngle, div)];
        idxVec{k} = [0 cumsum(iterVec{k})]+tempAngleVec(1)-1;
    end

    % generate colinear and transverse signals using interpolation
    x = gpuArray(-FOV_x/2:dx:FOV_x/2);
    z = gpuArray(-FOV_z/2:dz:FOV_z/2);
    colinearSignal = gpuArray.zeros(1, length(FFP_angle));
    transverseSignal = gpuArray.zeros(1, length(FFP_angle));
    horizontalSignal = gpuArray.zeros(size(colinearSignal));
    verticalSignal = gpuArray.zeros(size(colinearSignal));

    colinearIMG = gpuArray.zeros(zL, xL, div);
    transverseIMG = gpuArray.zeros(zL, xL, div);
    for j=1:length(numIters)       
        for l=1:numIters{j}
            colinearIMG(:,:,1:iterVec{j}(l)) = (h5read('./temp/PSF.h5','/colinearIMG', [1 1 idxVec{j}(l)+1], [zL xL iterVec{j}(l)]));
            transverseIMG(:,:,1:iterVec{j}(l)) = (h5read('./temp/PSF.h5','/transverseIMG', [1 1 idxVec{j}(l)+1], [zL xL iterVec{j}(l)]));
            try
                partialAngle = angleVec(idxVec{j}(l)+1:idxVec{j}(l+1));
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
            wait(gpudev); clear colinearIMG transverseIMG;
        end
        
    end

    % create relaxed colinear and transverse signals. Crop the signals 
    % appropriately.
    horizontalSignal_temp = conv([horizontalPrev horizontalSignal], r_t);
    horizontalPrev = horizontalSignal(end-2*Simparams.downsample-length(r_t)+1:end-2*Simparams.downsample);
    horizontalSignal = horizontalSignal_temp(length(horizontalPrev)+1:end-length(r_t)+1);

    verticalSignal_temp = conv([verticalPrev, verticalSignal], r_t);
    verticalPrev = verticalSignal(end-2*Simparams.downsample-length(r_t)+1:end-2*Simparams.downsample);
    verticalSignal = verticalSignal_temp(length(verticalPrev)+1:end-length(r_t)+1);
    

    SPIOparams.horizontalPrev{particleNo} = horizontalPrev;
    SPIOparams.verticalPrev{particleNo} = verticalPrev;

    % truncuate signals so that they have the correct sampling frequency (MPI frequency)
    truncSig = colinearSignal(:, 1:Simparams.downsample:end);
    signal.colinearSignal = gather(truncSig);
    truncSig = transverseSignal(:, 1:Simparams.downsample:end);
    signal.transverseSignal = gather(truncSig);
    truncSig = horizontalSignal(:, 1:Simparams.downsample:end);
    signal.horizontalSignal = gather(truncSig);
    truncSig = verticalSignal(:, 1:Simparams.downsample:end);
    signal.verticalSignal = gather(truncSig);
    
end
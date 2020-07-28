function generatePSFe(gpudev, MPIparams, SPIOparams, angleVec, particleNo, div)


    % create temporary data folder if it does not exist
    if ~exist('temp', 'dir')
           mkdir('temp')
    end
    % delete contents of temporary data folder
    tempList = dir('temp');
    if length(tempList) > 2
        delete('temp\*')
    end
    % divide data into chunks
    numAngle = length(angleVec);
    numIters = floor(numAngle/div);
    if ((numAngle/div) > numIters)
        numIters = numIters + 1;
    end
    if mod(numAngle, div) == 0
        addi = div;
    else
        addi = mod(numAngle, div);
    end
    iterVec = [div*ones(1, numIters-1), addi];
    idxVec = [0 cumsum(iterVec)];

    % compute PSFs
    G = [MPIparams.Gxx, MPIparams.Gyy, MPIparams.Gzz];
    FOV_x = MPIparams.FOV_x; 
    FOV_z = MPIparams.FOV_z;
    diameter = SPIOparams.diameter(particleNo);
    dx = SPIOparams.dx; 
    dz = SPIOparams.dz;
    
    x = gpuArray(-FOV_x/2:dx:FOV_x/2);
    z = gpuArray(-FOV_z/2:dz:FOV_z/2);
    xL = length(x); zL = length(z);
    % generate h5 files to store PSFs (this may be a big file depending on your parameters)
    h5create(['./temp/','PSF.h5'],'/colinearIMG',[zL xL numAngle]);
    h5create(['./temp/','PSF.h5'],'/transverseIMG',[zL xL numAngle]);
    h5create(['./temp/','PSF.h5'],'/angleVec',[1 length(angleVec)]);
    
    [X,Z] = meshgrid(x,z);
    
    mo = 4*pi*1e-7;
    Gxx = gpuArray(G(1)/mo);
    Gzz = gpuArray(G(3)/mo);

    d = diameter*1e-9;
    kb = 1.3806488e-23;
    T = 300;
    k = (0.1*pi*d^3)/(kb*T);
  
    e = gpuArray.zeros(2, xL*zL);
    e(1, :) = X(:);
    e(2, :) = Z(:);
    wait(gpudev); clear X Z;
    
    tempDistribution = SPIOparams.SPIOdistribution(:,:,particleNo);

    % pad the distribution with 0s so it has same size with the PSFs
    z_mpi = -MPIparams.FOV_z/2:dz:MPIparams.FOV_z/2;
    x_mpi = -MPIparams.FOV_x/2:dx:MPIparams.FOV_x/2;
    [xL_mpi, zL_mpi] = meshgrid(x_mpi, z_mpi);
    z_spio = -SPIOparams.image_FOV_z/2:dz:SPIOparams.image_FOV_z/2-dz;
    x_spio = -SPIOparams.image_FOV_x/2:dx:SPIOparams.image_FOV_x/2-dx;
    [xL_spio, zL_spio] = meshgrid(x_spio, z_spio);
    
    tempDistribution = interp2(xL_spio, zL_spio, tempDistribution,  xL_mpi, zL_mpi,'linear', 0);

    % multiplication in frequency domain is faster than convolution in
    % time domain
    imgF = gpuArray(fft2(tempDistribution, zL, xL));
    
    for l=1:numIters
        uniqueAngle_part = angleVec(idxVec(l)+1:idxVec(l+1));

        R = gpuArray.zeros(2,2, length(uniqueAngle_part));
        R(1,1, :) = cosd(uniqueAngle_part);
        R(1,2, :) = sind(uniqueAngle_part);
        R(2,1, :) = -sind(uniqueAngle_part);
        R(2,2, :) = cosd(uniqueAngle_part);

        C = pagefun(@mtimes,R,e); wait(gpudev); clear R;
        rotX = reshape(C(1, :), zL, xL, length(uniqueAngle_part));
        rotZ = reshape(C(2, :), zL, xL, length(uniqueAngle_part));
        wait(gpudev); clear C;
        Hxyzk = sqrt((Gxx*rotX).^2 + (Gzz*rotZ).^2)*k;
        idx = find(Hxyzk == 0);

        Lx = coth(Hxyzk) - 1./(Hxyzk) ;
        Lx = Lx./Hxyzk;
        Lx(idx) = 1/3; % Lx = Nenv, normal component from now on

        % transverse argument
        arg_trans = Gxx*Gzz*rotX.*rotZ./(Hxyzk/k).^2; arg_trans(idx) = 0;
        wait(gpudev); clear rotX;
        % collinear argument
        arg_colli = Gzz^2.*rotZ.^2./(Hxyzk/k).^2; arg_colli(idx) = 0;
        wait(gpudev); clear rotZ;

        Lx_der = 1./Hxyzk.^2 - 1./sinh(Hxyzk).^2;
        Lx_der(idx) = 1/3; % Lx_der = Tenv, tangential from now on
        wait(gpudev); clear Hxyzk;
        
        
        colinearPSF_f = fft2(Lx_der.*arg_colli + Lx.*(1-arg_colli)); wait(gpudev); 
        transversePSF_f = fft2((Lx_der - Lx).*arg_trans);
        wait(gpudev); clear Lx Lx_der arg_colli arg_trans;
        
        colIMG = real(ifft2(imgF.*colinearPSF_f)); wait(gpudev); clear colinearPSF_f;
        colIMG_s = fftshift(fftshift(colIMG, 1), 2); wait(gpudev); clear colIMG;
        h5write('./temp/PSF.h5', '/colinearIMG', gather(colIMG_s), [1 1 idxVec(l)+1], [zL xL iterVec(l)]); wait(gpudev); clear colIMG_s;
        
        tranIMG = real(ifft2(imgF.*transversePSF_f)); wait(gpudev); clear transversePSF_f;
        tranIMG_s = fftshift(fftshift(tranIMG, 1), 2); wait(gpudev); clear tranIMG;
        h5write('./temp/PSF.h5', '/transverseIMG', gather(tranIMG_s), [1 1 idxVec(l)+1], [zL xL iterVec(l)]); wait(gpudev); clear tranIMG_s;        
    end
    h5write('./temp/PSF.h5', '/angleVec', gather(angleVec), [1 1], [1 length(angleVec)]); wait(gpudev); clear angleVec; % write angle vector to file for ease of access later on

end
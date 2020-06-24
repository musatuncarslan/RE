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
    iterVec = [div*ones(1, numIters-1), mod(numAngle, div)];
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
    h5create(['./temp/','PSF.h5'],'/colinearPSF',[zL xL numAngle]);
    h5create(['./temp/','PSF.h5'],'/transversePSF',[zL xL numAngle]);
    
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
        Lx(idx) = 0;
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

        h5write('./temp/PSF.h5', '/colinearPSF', gather(Lx_der.*arg_colli + Lx.*(1-arg_colli)), [1 1 idxVec(l)+1], [zL xL iterVec(l)]);
        h5write('./temp/PSF.h5', '/transversePSF', gather((Lx_der - Lx).*arg_trans), [1 1 idxVec(l)+1], [zL xL iterVec(l)]);
    end

end
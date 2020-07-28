clear all
close all
clc

format long;

gpudev = gpuDevice(1); % get the GPU device

% parameters
Physicsparams = setPhysicsParams(); % physics parameters
MPIparams = setMPIParams(Physicsparams, [5, 0, 2.5]); % MPI machine parameters
SPIOparams = setSPIOParams(Physicsparams, 512, 2e-6); % SPIO parameters
[Simparams, MPIparams] = setSimulationParams(MPIparams, Physicsparams); % Simulation parameters

length(Simparams.simPeriods)

% preprocessing in order to find all the angles to calculate the necessary
% PSFs
maxIdx = 0;
for k=1:Simparams.divNum
    maxIdx = max(maxIdx, Simparams.div(k*2)-Simparams.div((k-1)*2+1)+1);
end
FFP_x = zeros(Simparams.divNum, maxIdx*Simparams.samplePerPeriod/Simparams.downsample+2);
FFP_z = zeros(size(FFP_x));
uniqueAngle_sim = [];
for k=1:Simparams.divNum
    idx = [Simparams.div((k-1)*2+1) Simparams.div(k*2)];
    simIdx = ((Simparams.simPeriods(idx(1))-1)*Simparams.samplePerPeriod+1:(Simparams.simPeriods(idx(2))*Simparams.samplePerPeriod+2*Simparams.downsample+1));
    t = gpuArray(simIdx/Simparams.fs_phsy);
    FFPparams = generateFFP(gpudev, t, MPIparams); 
    uniqueAngle_sim = [uniqueAngle_sim, FFPparams.FFP_uniqueAngle];
    FFP_x(k, 1:(idx(2)-idx(1)+1)*Simparams.samplePerPeriod/Simparams.downsample+2) = gather(FFPparams.FFP_x(1:Simparams.downsample:end-1));
    FFP_z(k, 1:(idx(2)-idx(1)+1)*Simparams.samplePerPeriod/Simparams.downsample+2) = gather(FFPparams.FFP_z(1:Simparams.downsample:end-1));
end
angleVec = unique(uniqueAngle_sim);

% plotting trajectory and SPIOs for check
figure; hold on;
x_axis = (-SPIOparams.image_FOV_x/2:SPIOparams.dx:SPIOparams.image_FOV_x/2-SPIOparams.dx);
z_axis = (-SPIOparams.image_FOV_z/2:SPIOparams.dz:SPIOparams.image_FOV_z/2-SPIOparams.dz);
for k=1:length(SPIOparams.diameter)
    surf(x_axis, z_axis, SPIOparams.SPIOdistribution(:,:,k)); shading interp
    hold on;
end
for k=1:Simparams.divNum
    idx = [Simparams.div((k-1)*2+1) Simparams.div(k*2)];
    plot(FFP_x(k, 1:(idx(2)-idx(1)+1)*Simparams.samplePerPeriod/Simparams.downsample+2), FFP_z(k, 1:(idx(2)-idx(1)+1)*Simparams.samplePerPeriod/Simparams.downsample+2), 'color', [0 0.4470 0.7410])
end
view(2)
xlabel('x-axis'); ylabel('z-axis')
xlim([-MPIparams.FOV_x/2 MPIparams.FOV_x/2]); ylim([-MPIparams.FOV_z/2 MPIparams.FOV_z/2])

wait(gpudev); clear FFPparams uniqueAngle_sim t;

% actual signal generation
maxIdx = 0;
for k=1:Simparams.divNum
    maxIdx = max(maxIdx, Simparams.div(k*2)-Simparams.div((k-1)*2+1)+1);
end
signal = zeros(Simparams.divNum, maxIdx*Simparams.samplePerPeriod/Simparams.downsample+2);
sigL = size(signal);
FFP_x = zeros(sigL);
FFP_z = zeros(sigL);
drive = zeros(sigL);
div = Simparams.div;
for particleNo = 1:length(SPIOparams.diameter)
    % efficient functions, uses HDD instead of RAM, should not give any
    % memory errors, uses HDF5 file structure
    tic
    generatePSFe(gpudev, MPIparams, SPIOparams, angleVec, particleNo, 40); 
    toc

    % generate the MPI signal
    for k=1:Simparams.divNum
        tic
        idx = [div((k-1)*2+1) div(k*2)]; % get the start and end indices of the necessary periods to be simulated, then the indices for each sample is generated
        simIdx = ((Simparams.simPeriods(idx(1))-1)*Simparams.samplePerPeriod+1:(Simparams.simPeriods(idx(2))*Simparams.samplePerPeriod+2*Simparams.downsample+1));
        t = gpuArray(simIdx/Simparams.fs_phsy); % divide the indices for each sample to sampling frequency to get the time vector
        FFPparams = generateFFP(gpudev, t, MPIparams); % simulate FFP
        wait(gpudev); clear t;
        [signals_sep, SPIOparams] = generateSe(gpudev, FFPparams, MPIparams, SPIOparams, Simparams, particleNo, 50); % simulate signal
        
        signal(k, 1:length(signals_sep.horizontalSignal)) = signal(k,1:length(signals_sep.horizontalSignal)) + gather(signals_sep.horizontalSignal); 
        wait(gpudev); clear t;
        
        if particleNo == 1
            FFP_x(k, 1:length(signals_sep.horizontalSignal)) = gather(FFPparams.FFP_x(1:Simparams.downsample:end-1));
            FFP_z(k, 1:length(signals_sep.horizontalSignal)) = gather(FFPparams.FFP_z(1:Simparams.downsample:end-1));
            drive(k, 1:length(signals_sep.horizontalSignal)) = gather(FFPparams.drive(1:Simparams.downsample:end-1));
        end
        wait(gpudev); clear signals_sep FFPparams;
        toc
    end
end

% delete signal file inside the temporary folder if it exists
tempList = dir('temp');
if length(tempList) > 3
    delete('temp\signal.h5')
end


h5create(['./temp/','signal.h5'],'/signal',sigL);
h5create(['./temp/','signal.h5'],'/FFP_x',sigL);
h5create(['./temp/','signal.h5'],'/FFP_z',sigL);
h5create(['./temp/','signal.h5'],'/drive',sigL);

h5write('./temp/signal.h5', '/signal', signal, [1 1], sigL);
h5write('./temp/signal.h5', '/FFP_x', FFP_x, [1 1], sigL);
h5write('./temp/signal.h5', '/FFP_z', FFP_z, [1 1], sigL);
h5write('./temp/signal.h5', '/drive', drive, [1 1], sigL);
% tic
% % regular PSF generation, uses RAM instead of HDD, may give memory errors
% % if angleVec is too long (i.e. angle resolution is high)
% [colinearPSF1,transversePSF1] = generatePSF(MPIparams, SPIOparams, FFPparams);
% [colinearIMG, transverseIMG] = generatePSFImages(SPIOparams, colinearPSF1, transversePSF1);
% toc

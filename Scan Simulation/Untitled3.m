clear all
close all
clc

format long;
format compact;

gpudev = gpuDevice(1); % get the GPU device

% parameters
Physicsparams = setPhysicsParams(); % physics parameters
MPIparams = setMPIParams(Physicsparams, [0.5, 0, 0.1]); % MPI machine parameters
SPIOparams = setSPIOParams(Physicsparams, 512, 2e-6); % SPIO parameters
[Simparams, MPIparams] = setSimulationParams(MPIparams, Physicsparams); % Simulation parameters


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
xlim([-0.05/2 0.05/2]); ylim([-0.06/2 0.06/2])

wait(gpudev); clear FFPparams uniqueAngle_sim t;

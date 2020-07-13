clear all
close all
clc

format long;
set(0,'FormatSpacing', 'compact')

gpudev = gpuDevice(1); % get the GPU device

addpath('./parameterFiles')
addpath('./signalGenerationFiles')

% parameters
Physicsparams = setPhysicsParams(); % physics parameters
MPIparams = setMPIParams(Physicsparams, 'complex_rastered', 0.10); % MPI machine parameters
SPIOparams = setSPIOParams(Physicsparams, 256, 2e-6); % SPIO parameters
[Simparams, MPIparams] = setSimulationParams(MPIparams, Physicsparams); % Simulation parameters



div = Simparams.div;
divNum = Simparams.divNum;

figure; plot(Simparams.simPeriods); 

figure;
%% partial trajectory 1
zSpeed = 1/(MPIparams.Gzz*(1/MPIparams.Rsz));
t = (Simparams.simPeriods-1)/MPIparams.f_drive;
p = 1/(MPIparams.numTrianglePeriods/MPIparams.time);
xPeriodsPosSim = MPIparams.FOV_x*(2*abs(2*(t/p - floor(t/p + 0.5)))-1)/2;
zPeriodsPosSim = t*zSpeed - MPIparams.FOV_z/2;
plot(xPeriodsPosSim, zPeriodsPosSim, 'd', 'color',[0.4660 0.6740 0.1880]); hold on

%% partial trajectory 2
for k=1:divNum
    tic
    idx = [div((k-1)*2+1) div(k*2)]; % get the start and end indices of the necessary periods to be simulated, then the indices for each sample is generated   
    simIdx = (Simparams.simPeriods(idx(1))-1)*Simparams.samplePerPeriod:(Simparams.simPeriods(idx(2))*Simparams.samplePerPeriod);
    t = gpuArray(simIdx/Simparams.fs_phsy); % divide the indices for each sample to sampling frequency to get the time vector
    zSpeed = 1/(MPIparams.Gzz*(1/MPIparams.Rsz));
    p = 1/(MPIparams.numTrianglePeriods/MPIparams.time);
    xPeriodsPosSim = MPIparams.FOV_x*(2*abs(2*(t/p - floor(t/p + 0.5)))-1)/2;
    zPeriodsPosSim = t*zSpeed - MPIparams.FOV_z/2;
    plot(xPeriodsPosSim, zPeriodsPosSim, '*', 'color', [0 0.4470 0.7410]); hold on
    toc
end



%% total trajectory
t = (0:MPIparams.zPeriods-1)/MPIparams.f_drive;
p = 1/(MPIparams.numTrianglePeriods/MPIparams.time);
xPeriodsPos = MPIparams.FOV_x*(2*abs(2*(t/p - floor(t/p + 0.5)))-1)/2;
zPeriodsPos = t*zSpeed - MPIparams.FOV_z/2;

plot(xPeriodsPos, zPeriodsPos, 'color', [0.6350 0.0780 0.1840])
xlim([-MPIparams.FOV_x/1.5 MPIparams.FOV_x/1.5]); ylim([-MPIparams.FOV_z/1.5 MPIparams.FOV_z/1.5]);
legend('Partial Trajectory', 'Complete Trajectory')
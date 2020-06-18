clear all
close all
clc

format long;

% parameters
Physicsparams = setPhysicsParams(); % physics parameters
MPIparams = setMPIParams(Physicsparams, 'complex_rastered', 0.1); % MPI machine parameters

SPIOparams = setSPIOParams(Physicsparams, 512, 2e-6); % SPIO parameters
Simparams = setSimulationParams(MPIparams, Physicsparams); % Simulation parameters





% simPeriodsx = findPeriodsX(MPIparams);
% simPeriodsz = findPeriodsZ(MPIparams);
% simPeriods = intersect(simPeriodsx,simPeriodsz,'stable'); % periods to simulate





% signal generation
time = MPIparams.time; % time to traverse whole FOV


numTrianglePeriods = MPIparams.numTrianglePeriods/MPIparams.time;
p = 1/numTrianglePeriods;
x_partial = []; % zeros(length(simPeriods)*samplePerPeriods, 1);
z_partial = []; % zeros(length(simPeriods)*samplePerPeriods, 1);
for k=Simparams.simPeriods
    tic
    t = (k*Simparams.samplePerPeriod:((k+1)*Simparams.samplePerPeriod-1))/Simparams.fs; 
    x_partial = [x_partial MPIparams.FOV_x*(2*abs(2*(t/p - floor(t/p + 0.5)))-1)/2]; % robot arm movement in x direction w.r.t. time
    z_partial = [z_partial MPIparams.FOV_z/time*t - MPIparams.FOV_z/2]; % robot arm movement in z direction w.r.t. time
    toc
end


t = (0:Simparams.fs*time-1)/Simparams.fs;
x = MPIparams.FOV_x*(2*abs(2*(t/p - floor(t/p + 0.5)))-1)/2;
z = MPIparams.FOV_z/time*t - MPIparams.FOV_z/2;

figure; plot(z_partial, x_partial, '*'); hold on;
plot(z, x, '.-');
xlabel('z-axis'); ylabel('x-axis')
legend('Partial Trajectory', 'Complete Trajectory')
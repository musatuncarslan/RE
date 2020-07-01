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

divNo = find(diff(Simparams.simPeriods) ~= 1); % index of division
if (isempty(divNo) ~= 1)
    divNo = [0 divNo length(Simparams.simPeriods)];
    divNum = length(divNo)-1; % number of divions
    divL = [];
    for k=1:divNum
        divL = [divL divNo(k)+1 divNo(k+1)] % length of each division
    end
else
    divNum = 1;
    divL = [1 length(Simparams.simPeriods)];
end
numTrianglePeriods = MPIparams.numTrianglePeriods/MPIparams.time;
p = 1/numTrianglePeriods;
x_partial = []; % zeros(length(simPeriods)*samplePerPeriods, 1);
z_partial = []; % zeros(length(simPeriods)*samplePerPeriods, 1);
for k=1:divNum
    idx = [divL((k-1)*2+1) divL(k*2)];
    simIdx = ((Simparams.simPeriods(idx(1))-1)*Simparams.samplePerPeriod+1:(Simparams.simPeriods(idx(2))*Simparams.samplePerPeriod+2*Simparams.downsample+1));
    t = gpuArray(simIdx/Simparams.fs_phsy);
    x_partial = [x_partial MPIparams.FOV_x*(2*abs(2*(t/p - floor(t/p + 0.5)))-1)/2]; % robot arm movement in x direction w.r.t. time
    z_partial = [z_partial MPIparams.FOV_z/time*t - MPIparams.FOV_z/2]; % robot arm movement in z direction w.r.t. time
end


t = (0:Simparams.fs_mpi*time-1)/Simparams.fs_mpi;
x = MPIparams.FOV_x*(2*abs(2*(t/p - floor(t/p + 0.5)))-1)/2;
z = MPIparams.FOV_z/time*t - MPIparams.FOV_z/2 + MPIparams.driveMag*cos(2*pi*MPIparams.f_drive*t);

figure;  
p1=plot3(x*100, z*100, zeros(1, length(z)), 'DisplayName', 'FFP Movement + Drive field'); hold on;
p2=plot3(x_partial*100, z_partial*100, zeros(1, length(z_partial)), 'DisplayName', 'FFP Movement');
xlim([0 0.02]); ylim([-1 1])
xlabel('x-axis (cm)'); ylabel('z-axis (cm)'); zlabel('y-axis (cm)')
legend('Location','northeast')
title({'Total FFP Movement under Drive Field (zoomed)', ['f_d = ' num2str(MPIparams.f_drive*1e-3) ' kHz']})
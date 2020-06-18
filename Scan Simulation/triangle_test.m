clear all
close all
clc

format long;

% parameters
Physicsparams = setPhysicsParams(); % physics parameters
MPIparams = setMPIParams(Physicsparams, 'complex_rastered', 0.1); % MPI machine parameters
SPIOparams = setSPIOParams(Physicsparams, 512, 2e-6); % SPIO parameters
Simparams = setSimulationParams(MPIparams, Physicsparams); % Simulation parameters



% % signal generation
% time = MPIparams.time; % time to traverse whole FOV
% numTrianglePeriods = MPIparams.numTrianglePeriods/MPIparams.time;
% p = 1/numTrianglePeriods;
% x_partial = []; % zeros(length(simPeriods)*samplePerPeriods, 1);
% z_partial = []; % zeros(length(simPeriods)*samplePerPeriods, 1);
% for k=Simparams.simPeriods
%     tic
%     t = (k*Simparams.samplePerPeriod:((k+1)*Simparams.samplePerPeriod-1))/Simparams.fs_phsy; 
%     x_partial = [x_partial MPIparams.FOV_x*(2*abs(2*(t/p - floor(t/p + 0.5)))-1)/2]; % robot arm movement in x direction w.r.t. time
%     z_partial = [z_partial MPIparams.FOV_z/time*t - MPIparams.FOV_z/2]; % robot arm movement in z direction w.r.t. time
%     toc
% end

x_sim = [];
z_sim = [];
speed_sim = [];
angle_sim = [];
uniqueAngle_sim = [];
count = 1;
sig_vert = [];
sig_hori = [];
for k=Simparams.simPeriods   

    tic
    t = gpuArray(((k-1)*Simparams.samplePerPeriod+1:(k*Simparams.samplePerPeriod+2*Simparams.downsample+1))/Simparams.fs_phsy);
    FFPparams = generateFFP(t, MPIparams, Simparams, [3, 0]);  
    x_sim = [x_sim, FFPparams.FFP_x(1:Simparams.downsample:end-1)];
    z_sim = [z_sim, FFPparams.FFP_z(1:Simparams.downsample:end-1)];
    speed_sim = [speed_sim, FFPparams.FFP_speed];
    angle_sim = [angle_sim, FFPparams.FFP_angle];
    uniqueAngle_sim = [uniqueAngle_sim, FFPparams.FFP_uniqueAngle];

%     in a single period, necessary angles will probably be gathered so
%     generating PSFs and corresponding images for a single period should
%     be enough
    if count == 1
        
        % calculate colinear and trasnverse PSF(s) for each unique angle
        [colinearPSF, transversePSF] = generatePSF(MPIparams, SPIOparams, FFPparams);
        
        
        % pad the empty parts of the image so that it has same size with PSF(s),
        % then compute colinear and transverse images for each unique angle
        
        [colinearIMG, transverseIMG] = generatePSFImages(SPIOparams, colinearPSF, transversePSF);
        clearvars colinearPSF transversePSF;
    end
    % generate MPI signals
    [signals_sep, SPIOparams] = generateSignals(colinearIMG, transverseIMG, FFPparams, MPIparams, SPIOparams, Simparams);
    
    signals = struct;
    signals.horizontalSignal = zeros(1, size(signals_sep, 2));
    signals.verticalSignal = zeros(1, size(signals_sep, 2));
    for l=1:length(SPIOparams.diameter)
        signals.horizontalSignal = signals.horizontalSignal + real(signals_sep.horizontalSignal(l, :));
        signals.verticalSignal = signals.verticalSignal + real(signals_sep.verticalSignal(l, :));
    end
    
    sig_vert = [sig_vert, signals.verticalSignal];
    sig_hori = [sig_hori, signals.horizontalSignal];
    
    if k == Simparams.simPeriods(end/2)
        plot(signals.horizontalSignal(1:100)); hold on; plot(-flip(signals.horizontalSignal(102:end)))
        xlabel('Sample #'); ylabel('s(t)'); legend('s_{pos}(t)', '-flip(s_{pos}(t)')
    end
    
%     plot(t(1:Simparams.downsample:end-1), signals.horizontalSignal, 'blue'); hold on;
    
    count = count + 1;
    toc
end

figure; scatter3(x_sim, z_sim, sig_hori, 3, sig_hori)
title('Impulse SPIO @ center')
xlabel('z-axis (m)'); ylabel('x-axis (m)'); zlabel('s(t)')


p = 1/(MPIparams.numTrianglePeriods/MPIparams.time);
t = (0:Simparams.fs_mpi*MPIparams.time-1)/Simparams.fs_mpi;
x = MPIparams.FOV_x*(2*abs(2*(t/p - floor(t/p + 0.5)))-1)/2;
z = MPIparams.FOV_z/MPIparams.time*t - MPIparams.FOV_z/2;

figure; hold on;
plot(z, x, '.-');
plot(z_sim, x_sim)
xlabel('z-axis'); ylabel('x-axis')
legend('Complete Trajectory', 'Partial Trajectory')
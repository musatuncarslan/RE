clear all
close all
clc

format long;

addpath('./parameterFiles')
addpath('./signalGenerationFiles')

% parameters
Physicsparams = setPhysicsParams(); % physics parameters
MPIparams = setMPIParams(Physicsparams, 20); % MPI machine parameters
% SPIOparams = setSPIOParams(Physicsparams, 256, 2e-6); % SPIO parameters
Simparams = setSimulationParams(MPIparams, Physicsparams); % Simulation parameters

slewRate = [0.0001 0.25:0.25:20]; % T/s
Hp=MPIparams.Bp/Physicsparams.mu0; % magnetization moment
G=MPIparams.Gzz/Physicsparams.mu0; % gradient
MPIparams.driveMag=Hp/G; % extent of the drive field

R = slewRate/MPIparams.Gzz;


f_d = MPIparams.f_drive; 
fs = 100e6;
t = (0:fs*5/f_d-1)/fs; % time axis for error calculation



for l=1:length(R)
    MPIparams.slewRate = slewRate(l);
    [del_t(l), a_t(l)] = findDistortion(MPIparams);
    del_t(l) = del_t(l) + 1/2/MPIparams.f_drive; 
    err(l) = abs(MPIparams.driveMag*cos(2*pi*f_d/4/f_d) - (MPIparams.driveMag*cos(2*pi*f_d*(1/(4*f_d) + del_t(l))) + R(l)*del_t(l)));

end
figure; 
subplot(2,2,1); plot(slewRate, (del_t-1/(2*f_d))*1e6); title('\Delta t w.r.t. Slew Rate'); xlabel('R_s (T/s)'); ylabel('\Delta t (\mu s)'); set(gca,'FontSize',12)
subplot(2,2,2); plot(slewRate, err*10^6); title('Error of Analytic solution'); xlabel('R_s (T/s)'); ylabel({'Absolute Error', 'in Position (\mu m)'}); ylim([0 1e-3]); set(gca,'FontSize',12)

subplot(2,2,[3, 4])
idx = length(R);
plot(t, (MPIparams.driveMag*cos(2*pi*f_d*t) + R(idx)*t)*1e3); hold on;
plot(1/4/f_d, (MPIparams.driveMag*cos(2*pi*f_d/4/f_d) + R(idx)/4/f_d)*1e3, '*');
plot(1/4/f_d+del_t(idx), (MPIparams.driveMag*cos(2*pi*f_d*(1/4/f_d+del_t(idx))) + R(idx)*(1/4/f_d+del_t(idx)))*1e3, '*');
set(gca,'FontSize',12)
title(['FFP Position w.r.t. time for R_s = ' num2str(slewRate(idx)) ' T/s'])
ylabel('z-axis position (mm)'); xlabel('t (s)'); legend('FFP movement', 'Position @ max speed', 'Analytic sol. to pos. @ max speed')

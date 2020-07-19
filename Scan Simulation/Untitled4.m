clear all
close all
clc

format long
format compact

% parameters
Physicsparams = setPhysicsParams(); % physics parameters
MPIparams = setMPIParams(Physicsparams, 'complex_rastered', 0.1); % MPI machine parameters

Rsz = 0.1*ones(1, 10);
Rsx = [0.66, 1, 1.33, 1.66, 2, 2.33, 4, 6.66, 13.33, 26.66];

meanVal = [3.75, 3.74, 3.76, 3.74, 3.75, 3.74, 3.75, 3.80, 3.88, 4.06];
stdVal = [0.036,0.082,0.098,0.119,0.145,0.167,0.280,0.467,0.974,1.702];
nRMSE = [0.013,0.022,0.026,0.033,0.039,0.045,0.075,0.123,0.250,0.419];

Rsx_interp = [0.66:0.01:20];
nRMSE_interp = interp1(Rsx, nRMSE, Rsx_interp, 'spline');


xd_interp = Rsx_interp/4.8/10^4;
figure; plot(xd_interp./MPIparams.driveMag, nRMSE_interp); hold on;
title({'nRMSE vs. Normalized Displacement Distance in x-axis', '$R_{s,z} = 0.1$ T/s'}, 'interpreter', 'latex')
xlabel('$x_d/pFOV$', 'interpreter', 'latex'); ylabel('nRMSE')

figure; plot(Rsx_interp, nRMSE_interp); hold on;
title({'nRMSE vs. Normalized Displacement Distance in x-axis', '$R_{s,z} = 0.1$ T/s'}, 'interpreter', 'latex')
xlabel('$R_{s,x}$ T/s', 'interpreter', 'latex'); ylabel('nRMSE')
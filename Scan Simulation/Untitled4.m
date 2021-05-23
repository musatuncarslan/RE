clear all
close all
clc

format long
format compact

% parameters
Physicsparams = setPhysicsParams(); % physics parameters
MPIparams = setMPIParams(Physicsparams, [0 0 3]); % MPI machine parameters

fileID = fopen('last_data.txt');
C = textscan(fileID,'%f %f %f %f %f %f %*[^\n]', 'Delimiter', ' ');
fclose(fileID);

datL = 20*10;
skipIdx = []; % [6*12+1:7*12, 8*12+1:9*12];

Rsx = (C{1}(1:datL)); Rsx(skipIdx) = [];
Rsz = C{2}(1:datL); Rsz(skipIdx) = [];
meanVal = C{3}(1:datL); meanVal(skipIdx) = [];
stdVal = C{4}(1:datL); stdVal(skipIdx) = [];
nRMSE = C{5}(1:datL); nRMSE(skipIdx) = [];

Rsx_interp = 0.1:0.1:10;
[Rsx_interp, Rsz_interp] = meshgrid(Rsx_interp, linspace(0.1, max(Rsz), length(Rsx_interp)));


F = scatteredInterpolant(Rsx,Rsz,100*(1-(meanVal/4)));
F.Method = 'natural';
mean_interp = F(Rsx_interp,Rsz_interp);
% nRMSE_interp = interp2(Rsx, Rsz, nRMSE, Rsx_interp, Rsz_interp, 'spline');


figure; surf(Rsx_interp, Rsz_interp, mean_interp);  shading interp; hold on;
title({'Estimation Error (\%) vs. Slew Rate', 'with Time Shift Correction'}, 'interpreter', 'latex')
xlabel('$R_{s,x}$', 'interpreter', 'latex', 'fontsize', 14); ylabel('$R_{s,z}$', 'interpreter', 'latex', 'fontsize', 14)
zlabel('nRMSE')
axis tight; view(2)
colormap(parula(datL*0.1))
hcb = colorbar;
set(get(hcb,'Title'),'String','\%','interpreter', 'latex')

% figure; plot(Rsx_interp, nRMSE_interp); hold on;
% title({'nRMSE vs. Normalized Displacement Distance in x-axis', '$R_{s,z} = 0.1$ T/s'}, 'interpreter', 'latex')
% xlabel('$R_{s,x}$ T/s', 'interpreter', 'latex'); ylabel('nRMSE')
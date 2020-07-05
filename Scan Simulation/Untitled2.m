clear all
close all
clc

format long;
set(0,'FormatSpacing', 'compact')

addpath('./parameterFiles')
addpath('./signalGenerationFiles')

% parameters
Physicsparams = setPhysicsParams(); % physics parameters
MPIparams = setMPIParams(Physicsparams, 'complex_rastered', 0.1); % MPI machine parameters
% SPIOparams = setSPIOParams(Physicsparams, 256, 2e-6); % SPIO parameters
% Simparams = setSimulationParams(MPIparams, Physicsparams); % Simulation parameters

MPIparams.f_drive = 10.800e3;
zspeed = 1/(MPIparams.Gzz*(1/MPIparams.slewRate));
x = floor(MPIparams.FOV_z*MPIparams.f_drive/zspeed) % integer number of periods in the whole fov (this changes the FOV_z distance slightly)
FOV_z_a = (x)/MPIparams.f_drive*zspeed % adjusted z-FOV
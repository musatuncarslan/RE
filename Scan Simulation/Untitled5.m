clear all
close all
clc

format long;
format compact;

gpudev = gpuDevice(1); % get the GPU device

rs_z = [0.1, 0.5, round(linspace(1, 20, 20), 2)];
rs_x = 0.0001; % [0.1, 0.5, round(linspace(1, 10, 10), 2)];

count = 1;
for rs_z_idx = rs_z
   for rs_x_idx = rs_x
        % parameters
        Physicsparams = setPhysicsParams(); % physics parameters
        MPIparams = setMPIParams(Physicsparams, [rs_x_idx, 0, rs_z_idx]); % MPI machine parameters
        SPIOparams = setSPIOParams(Physicsparams, 512, 2e-6); % SPIO parameters
%         [Simparams, MPIparams] = setSimulationParams(MPIparams, Physicsparams); % Simulation parameters

        % find estimation parameters for phase and amplitude correction
        estimationParams = struct;
        [del_t(count, 1), amp_t(count, 1)] = findDistortion(MPIparams); 
        count = count + 1;
   end
end


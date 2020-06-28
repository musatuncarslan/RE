clear all
close all
clc

format long;

gpudev = gpuDevice(1); % get the GPU device

% parameters
Physicsparams = setPhysicsParams(); % physics parameters
MPIparams = setMPIParams(Physicsparams, 'complex_rastered', 0.1); % MPI machine parameters
SPIOparams = setSPIOParams(Physicsparams, 512, 2e-6); % SPIO parameters
Simparams = setSimulationParams(MPIparams, Physicsparams); % Simulation parameters
particleNo = 1;





% load signal_08-June-2020_14-30-33;

info = h5info('./temp/signal.h5');

sigL = info.Datasets(1).Dataspace.Size(2);
horizontalSignal_mpi_mat = h5read('./temp/signal.h5', '/signal', [1 1], [1 sigL]);
FFP_z = h5read('./temp/signal.h5', '/FFP_z', [1 1], [1 sigL]);
FFP_x = h5read('./temp/signal.h5', '/FFP_x', [1 1], [1 sigL]);
FFP_drive = h5read('./temp/signal.h5', '/drive', [1 1], [1 sigL]);

% arrange matrix data into vector
numIters = 1;
numSamplesPerIter = sigL - 2;
horizontalSignal = zeros(1, numIters*numSamplesPerIter+2);
FFP_z_fixed = zeros(1, numIters*numSamplesPerIter+2);
FFP_drive_fixed = zeros(1, numIters*numSamplesPerIter+2);
FFP_speed_fixed = zeros(1, numIters*numSamplesPerIter+2);
for k=1:numIters
    if k ~= numIters
        horizontalSignal(1, numSamplesPerIter*(k-1)+1:numSamplesPerIter*k) = horizontalSignal_mpi_mat(k, 1:end-2);
        FFP_z_fixed(1, numSamplesPerIter*(k-1)+1:numSamplesPerIter*k) = FFP_z(k, 1:end-2);
        FFP_x_fixed(1, numSamplesPerIter*(k-1)+1:numSamplesPerIter*k) = FFP_x(k, 1:end-2);
        FFP_drive_fixed(1, numSamplesPerIter*(k-1)+1:numSamplesPerIter*k) = FFP_drive(k, 1:end-2);
    else 
        horizontalSignal(1, numSamplesPerIter*(k-1)+1:numSamplesPerIter*k+2) = horizontalSignal_mpi_mat(k, :);
        FFP_z_fixed(1, numSamplesPerIter*(k-1)+1:numSamplesPerIter*k+2) = FFP_z(k, :);
        FFP_x_fixed(1, numSamplesPerIter*(k-1)+1:numSamplesPerIter*k+2) = FFP_x(k, :);
        FFP_drive_fixed(1, numSamplesPerIter*(k-1)+1:numSamplesPerIter*k) = FFP_drive(k, 1:end-2);
    end
end

L = sigL;
f = (0:L-1)*(MPIparams.fs)/L;
f_step = MPIparams.fs/L;
Bw_idx = round(10000/f_step/2);

% extract odd harmonics
idx = [];
for k=3
    [~, idx_1] = min(abs(f-MPIparams.f_drive*k));
    [~, idx_2] = min(abs(f-(MPIparams.fs-MPIparams.f_drive*k)));
    idx_1 = idx_1-Bw_idx:idx_1+Bw_idx;
    idx_2 = idx_2-Bw_idx:idx_2+Bw_idx;
    idx = [idx idx_1 idx_2];
end


% find estimation parameters for phase and amplitude correction
estimationParams = struct;
[estimationParams.del_t, estimationParams.amp_t] = findDistortion(MPIparams); 

% filter the 1st harmonic and extract odd harmonics upto 7th. also upsample
% the data if necessary
[filteredSignal, MPIparams] = extractData(horizontalSignal, MPIparams, 2e6, 1000, 11);

snr = inf;
[signal, ~] = awgnInterference(filteredSignal, MPIparams, snr, inf);




% start processing
start_idx = 1; % start_idx+val_idx-3; % corrected start idx "-2" is magical dont ask
pfov = (sigL-2)/(Simparams.samplePerPeriod/Simparams.downsample); % total number of pfovs
p_start = 0; % floor((start_idx)/1000) get back to the start of the data from the starting sample
p_end = pfov-p_start;

data_start_idx = start_idx; % 5309152, 6912443+3700234+242-137 % 10368663+813+10-43; % 16922286-100+569+500-5+194
data_idx = data_start_idx-(p_start)*1e3:data_start_idx+(p_end)*MPIparams.fs/MPIparams.f_drive+1;
partial_signal_interp = signal(data_idx);

numIters = pfov/4; 
numPeriod = pfov; % (length(partial_signal_interp)-2)/(MPIparams.fs/MPIparams.f_drive); % number of periods on a single line

estimationParams.interp_coeff = 100;
numPeriodsPerIter = numPeriod/numIters;
estimationParams.numSamplesPerIter = 1*numPeriodsPerIter*(MPIparams.fs/MPIparams.f_drive)+2;
estimationParams.numSample = (MPIparams.fs/MPIparams.f_drive)+2; % + 2 is for interpolation reasons
estimationParams.numSampleInterpolated = estimationParams.numSamplesPerIter*estimationParams.interp_coeff;

tau_est_linear = zeros(1, numIters);
tau_est_frequency = zeros(1, numIters);
tau_lin_weighted = zeros(1, numIters);
count = 1;
for pfovIdx=1:numIters
    tic
    % sliding window
%     sig_idx = (1:(estimationParams.numSamplesPerIter-2)+2)+(estimationParams.numSample-2)*(pfovIdx-1); 
    % non-sliding window    
    sig_idx = (1:(estimationParams.numSamplesPerIter-2)+2)+(estimationParams.numSamplesPerIter-2)*(pfovIdx-1);
    FFP_z_part(count) = FFP_z_fixed(sig_idx(end/2)) - FFP_drive_fixed(sig_idx(end/2));
    FFP_x_part(count) = FFP_x_fixed(sig_idx(end/2));
    
    [tau_est_frequency(count), tau_est_linear(count), tau_lin_weighted(count)] = estimateTau_func(MPIparams, estimationParams, partial_signal_interp(sig_idx), 5);
%     err(count) = calculateError(MPIparams, estimationParams, partial_signal_interp(sig_idx), tau_est_linear(count));
    toc
    count = count + 1;
end

figure; scatter3(FFP_x_part, FFP_z_part, tau_est_linear*10^6,  3, tau_est_linear*10^6)
xlabel('x-axis'); ylabel('z-axis'); zlabel('\tau (us)')


% parameters
Physicsparams = setPhysicsParams(); % physics parameters
MPIparams = setMPIParams(Physicsparams, 'complex_rastered', 0.1); % MPI machine parameters
Physicsparams.fs = MPIparams.fs;
SPIOparams = setSPIOParams(Physicsparams, 512, 2e-6); % SPIO parameters
Simparams = setSimulationParams(MPIparams, Physicsparams); % Simulation parameters


figure;
x_axis = (-SPIOparams.image_FOV_x/2:SPIOparams.dx:SPIOparams.image_FOV_x/2-SPIOparams.dx);
z_axis = (-SPIOparams.image_FOV_z/2:SPIOparams.dz:SPIOparams.image_FOV_z/2-SPIOparams.dz);
for k=1:length(SPIOparams.diameter)
    surf(x_axis, z_axis, SPIOparams.SPIOdistribution(:,:,k)); shading interp
    hold on;
end
for k=1   
    t = gpuArray(((Simparams.simPeriods(1)-1)*Simparams.samplePerPeriod+1:(Simparams.simPeriods(end)*Simparams.samplePerPeriod+2*Simparams.downsample+1))/Simparams.fs_phsy);
    FFPparams = generateFFP(t, MPIparams, Simparams, [3, 0]); 
end
plot(FFPparams.FFP_x(1:Simparams.downsample:end-1), FFPparams.FFP_z(1:Simparams.downsample:end-1))
view(2)
xlabel('x-axis'); ylabel('z-axis')

F = griddedInterpolant({FFP_x_part', FFP_z_part'}, tau_est_linear');
F.Method = 'natural';
[xq, zq] = meshgrid(x_axis, z_axis);
Vq = F(xq, zq);


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


info = h5info('./temp/signal.h5');

sigSize = info.Datasets(1).Dataspace.Size;
horizontalSignal_mpi_mat = h5read('./temp/signal.h5', '/signal', [1 1], sigSize);
FFP_z = h5read('./temp/signal.h5', '/FFP_z', [1 1], sigSize);
FFP_x = h5read('./temp/signal.h5', '/FFP_x', [1 1], sigSize);
FFP_drive = h5read('./temp/signal.h5', '/drive', [1 1], sigSize);

% arrange matrix data into vector
numIters = sigSize(1);
numSamplesPerIter = sigSize(2) - 2;
horizontalSignal = zeros(1, numIters*numSamplesPerIter+2);
FFP_z_fixed = zeros(1, numIters*numSamplesPerIter+2);
FFP_drive_fixed = zeros(1, numIters*numSamplesPerIter+2);
FFP_x_fixed = zeros(1, numIters*numSamplesPerIter+2);
divL = Simparams.divL;
L = 0;
for k=1:numIters
    idx = [1 divL(k*2)-divL((k-1)*2+1)+1];
    if k ~= numIters
        vecIdx = ((idx(1)-1)*Simparams.samplePerPeriod/Simparams.downsample+1:idx(2)*Simparams.samplePerPeriod/Simparams.downsample)+L;
        horizontalSignal(1, vecIdx) = horizontalSignal_mpi_mat(k, 1:length(vecIdx));
        FFP_z_fixed(1, vecIdx) = FFP_z(k, 1:length(vecIdx));
        FFP_x_fixed(1, vecIdx) = FFP_x(k, 1:length(vecIdx));
        FFP_drive_fixed(1, vecIdx) = FFP_drive(k, 1:length(vecIdx));
    else 
        vecIdx = ((idx(1)-1)*Simparams.samplePerPeriod/Simparams.downsample+1:idx(2)*Simparams.samplePerPeriod/Simparams.downsample)+2+L;
        horizontalSignal(1, vecIdx) = horizontalSignal_mpi_mat(k, 1:length(vecIdx));
        FFP_z_fixed(1, vecIdx) = FFP_z(k, 1:length(vecIdx));
        FFP_x_fixed(1, vecIdx) = FFP_x(k, 1:length(vecIdx));
        FFP_drive_fixed(1, vecIdx) = FFP_drive(k, 1:length(vecIdx));
    end
    L = L+length(vecIdx);
end

figure; scatter3(FFP_x_fixed*100, FFP_z_fixed*100, horizontalSignal,  3, horizontalSignal)
xlabel('x-axis (cm)'); ylabel('z-axis (cm)'); zlabel('s(t)'); title('MPI Signal')


L = length(horizontalSignal);
f = (0:L-1)*(MPIparams.fs)/L;
f_step = MPIparams.fs/L;
Bw_idx = round(10000/f_step/2);

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
pfov = length(Simparams.simPeriods); % total number of pfovs
p_start = 0; % floor((start_idx)/1000) get back to the start of the data from the starting sample
p_end = pfov-p_start;

data_start_idx = start_idx; % 5309152, 6912443+3700234+242-137 % 10368663+813+10-43; % 16922286-100+569+500-5+194
data_idx = data_start_idx-(p_start)*1e3:data_start_idx+(p_end)*MPIparams.fs/MPIparams.f_drive+1;
partial_signal_interp = signal(data_idx);

numIters = pfov/8; 
numPeriod = pfov; % (length(partial_signal_interp)-2)/(MPIparams.fs/MPIparams.f_drive); % number of periods on a single line

estimationParams.interp_coeff = 100;
numPeriodsPerIter = numPeriod/numIters;
estimationParams.numSamplesPerIter = 1*numPeriodsPerIter*(MPIparams.fs/MPIparams.f_drive)+2;
estimationParams.numSample = (MPIparams.fs/MPIparams.f_drive)+2; % + 2 is for interpolation reasons
estimationParams.numSampleInterpolated = estimationParams.numSamplesPerIter*estimationParams.interp_coeff;

tau_est_linear = zeros(1, numIters);
tau_est_frequency = zeros(1, 692928);
tau_lin_weighted = zeros(1, numIters);
count = 1;
for pfovIdx=1:numIters
    tic
    % sliding window
%     sig_idx = (1:(estimationParams.numSamplesPerIter-2)+2)+(estimationParams.numSample-2)*(pfovIdx-1); 
    % non-sliding window    
    sig_idx = ((1:(estimationParams.numSamplesPerIter-2)+2)+(estimationParams.numSamplesPerIter-2)*(pfovIdx-1));
    FFP_z_part(count) = FFP_z_fixed(sig_idx(end/2)) - FFP_drive_fixed(sig_idx(end/2));
    FFP_x_part(count) = FFP_x_fixed(sig_idx(end/2));
    
    [tau_est_frequency(count), tau_est_linear(sig_idx), tau_lin_weighted(count)] = estimateTau_func(MPIparams, estimationParams, partial_signal_interp(sig_idx), 5);
%     err(count) = calculateError(MPIparams, estimationParams, partial_signal_interp(sig_idx), tau_est_linear(count));
    toc
    count = count + 1;
end

figure; scatter3(FFP_x_fixed, FFP_z_fixed, tau_est_linear*10^6,  3, tau_est_linear*10^6)
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
distribution = zeros(size(SPIOparams.SPIOdistribution(:,:,1)));
for k=1:length(SPIOparams.diameter)
    distribution = distribution + SPIOparams.SPIOdistribution(:,:,k);
end
surf(x_axis*100, z_axis*100, distribution); shading interp
hold on;
% preprocessing in order to find all the angles to calculate the necessary
% PSFs
maxIdx = 0;
for k=1:Simparams.divNum
    maxIdx = max(maxIdx, Simparams.divL(k*2)-Simparams.divL((k-1)*2+1)+1);
end
FFP_x = zeros(Simparams.divNum, maxIdx*Simparams.samplePerPeriod/Simparams.downsample+2);
FFP_z = zeros(size(FFP_x));
uniqueAngle_sim = [];
for k=1:Simparams.divNum
    idx = [Simparams.divL((k-1)*2+1) Simparams.divL(k*2)];
    simIdx = ((Simparams.simPeriods(idx(1))-1)*Simparams.samplePerPeriod+1:(Simparams.simPeriods(idx(2))*Simparams.samplePerPeriod+2*Simparams.downsample+1));
    t = gpuArray(simIdx/Simparams.fs_phsy);
    FFPparams = generateFFP(gpudev, t, MPIparams, Simparams, [3, 0]); 
    uniqueAngle_sim = [uniqueAngle_sim, FFPparams.FFP_uniqueAngle];
    plot(FFPparams.FFP_x(1:Simparams.downsample:end-1)*100, FFPparams.FFP_z(1:Simparams.downsample:end-1)*100, 'color', [0 0.4470 0.7410])
end
view(2)
xlabel('x-axis (cm)'); ylabel('z-axis (cm)'); title('Impulsive SPIO Distribution'); legend('SPIO distribution', 'Total FFP Trajectory')

F = scatteredInterpolant(FFP_x_fixed', FFP_z_fixed', tau_est_linear'*10^6);
F.Method = 'natural';
[xq, zq] = meshgrid(x_axis, z_axis);
Vq = F(xq, zq);

figure; surf(x_axis*100, z_axis*100, Vq.*distribution)
shading interp
xlabel('x-axis (cm)'); ylabel('z-axis (cm)'); zlabel('\tau (\mu s)')
axis tight; title(['Estimation Surface, f_d = ' num2str(MPIparams.f_drive*1e-3) ' kHz']); 
xlim([-1 1]); ylim([-1 1])
colorbar; 
% view(2);

meanTau = [];
stdTau = [];
for k=1:length(SPIOparams.diameter)
    meanTau = [meanTau mean(Vq(SPIOparams.SPIOdistribution(:,:,k)~=0))];
    stdTau = [stdTau std(Vq(SPIOparams.SPIOdistribution(:,:,k)~=0))];
end
meanTau
stdTau
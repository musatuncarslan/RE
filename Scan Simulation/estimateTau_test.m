format long;
format compact;

gpudev = gpuDevice(1); % get the GPU device

info = h5info('./temp/signal.h5');

sigSize = info.Datasets(1).Dataspace.Size;
horizontalSignal_mpi_mat = h5read('./temp/signal.h5', '/signal', [1 1], sigSize);
FFP_z = h5read('./temp/signal.h5', '/FFP_z', [1 1], sigSize);
FFP_x = h5read('./temp/signal.h5', '/FFP_x', [1 1], sigSize);
FFP_drive = h5read('./temp/signal.h5', '/drive', [1 1], sigSize);

% % arrange matrix data into vector
horizontalSignal = horizontalSignal_mpi_mat;
FFP_z_fixed = FFP_z;
FFP_x_fixed = FFP_x;
FFP_drive_fixed = FFP_drive;

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
interp_coeff = 1;
[filteredSignal, MPIparams] = extractData(horizontalSignal, MPIparams, MPIparams.fs*interp_coeff, 1000, 11);

snr = inf;
[signal, ~] = awgnInterference(filteredSignal, MPIparams, snr, inf);

% start processing
start_idx = 1; % start_idx+val_idx-3; % corrected start idx "-2" is magical dont ask
pfov = (length(filteredSignal)-2)/(Simparams.samplePerPeriod/Simparams.downsample*interp_coeff); % length(Simparams.simPeriods); % total number of pfovs
data_idx = 1:pfov*MPIparams.fs/MPIparams.f_drive+2;
partial_signal_interp = signal(data_idx);

numPeriodsPerIter = 4;
numIters = pfov/numPeriodsPerIter; 
numPeriod = pfov; % (length(partial_signal_interp)-2)/(MPIparams.fs/MPIparams.f_drive); % number of periods on a single line

estimationParams.interp_coeff = 1;
estimationParams.numSamplesPerIter = numPeriodsPerIter*(Simparams.fs_mpi/MPIparams.f_drive)*interp_coeff+2;
estimationParams.numSample = (MPIparams.fs/MPIparams.f_drive)+2; % + 2 is for interpolation reasons
estimationParams.numSampleInterpolated = estimationParams.numSamplesPerIter*estimationParams.interp_coeff;

tau_est_linear = zeros(1, numIters);
tau_est_frequency = zeros(1, numIters);
tau_lin_weighted = zeros(1, numIters);
count = 1;
for pfovIdx=round(numIters/2)
    % sliding window
%     sig_idx = (1:(estimationParams.numSamplesPerIter-2)+2)+(estimationParams.numSample-2)*(pfovIdx-1); 
    % non-sliding window    
    sig_idx = ((1:(estimationParams.numSamplesPerIter-2))+(estimationParams.numSamplesPerIter-2)*(pfovIdx-1));
    
    [tau_est_frequency(count), tau_est_linear(pfovIdx), tau_lin_weighted(sig_idx)] = estimateTau_func_new(MPIparams, estimationParams, filteredSignal(sig_idx), 5);
    count = count + 1;
end



% figure; scatter3(FFP_x_fixed, FFP_z_fixed, tau_est_linear*10^6,  3, tau_est_linear*10^6)
% xlabel('x-axis'); ylabel('z-axis'); zlabel('\tau (us)')
% 
% figure;
% x_axis = (-SPIOparams.image_FOV_x/2:SPIOparams.dx:SPIOparams.image_FOV_x/2-SPIOparams.dx);
% z_axis = (-SPIOparams.image_FOV_z/2:SPIOparams.dz:SPIOparams.image_FOV_z/2-SPIOparams.dz);
% distribution = zeros(size(SPIOparams.SPIOdistribution(:,:,1)));
% for k=1:length(SPIOparams.diameter)
%     distribution = distribution + SPIOparams.SPIOdistribution(:,:,k);
% end
% surf(x_axis*100, z_axis*100, distribution); shading interp
% title('Distribution')
% hold on;
% 
% F = scatteredInterpolant(FFP_x_fixed', FFP_z_fixed', tau_est_linear'*10^6);
% F.Method = 'natural';
% [xq, zq] = meshgrid(x_axis, z_axis);
% Vq = F(xq, zq);
% 
% figure; surf(x_axis*100, z_axis*100, Vq.*distribution)
% shading interp
% xlabel('x-axis (cm)'); ylabel('z-axis (cm)'); zlabel('\tau (\mu s)')
% axis tight; title(['Estimation Surface, f_d = ' num2str(MPIparams.f_drive*1e-3) ' kHz']); 
% colorbar; 
% view(2);
% 
% meanTau = [];
% stdTau = [];
% for k=1:length(SPIOparams.diameter)
%     meanTau = [meanTau mean(Vq(SPIOparams.SPIOdistribution(:,:,k)~=0))];
%     stdTau = [stdTau std(Vq(SPIOparams.SPIOdistribution(:,:,k)~=0))];
%     RMSE = sqrt(sum((Vq(SPIOparams.SPIOdistribution(:,:,k)~=0)-meanTau(k)).^2)/numel(find(SPIOparams.SPIOdistribution(:,:,k)~=0)));
% end
MPIparams.Rs(1);
MPIparams.Rs(3);
meanTau = tau_est_linear*10^6;
round(meanTau, 3);
stdTau = 0;
round(stdTau, 3);
RMSE = 0;
nRMSE = round(RMSE/meanTau, 3);
% nRMSE = round(RMSE/(max(Vq(SPIOparams.SPIOdistribution(:,:,k)~=0))-min(Vq(SPIOparams.SPIOdistribution(:,:,k)~=0))), 3);

fprintf('%1.2f %12.2f %12.3f %12.3f %12.3f\n', MPIparams.Rs(1), MPIparams.Rs(3), ...
    round(meanTau, 5), round(stdTau, 5), round(RMSE/meanTau, 5))

% 
% 

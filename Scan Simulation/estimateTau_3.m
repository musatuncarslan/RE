clear all
close all
clc

format long;

load signal_08-June-2020_14-30-33;



% arrange matrix data into vector
numIters = size(FFP_z, 1);
numSamplesPerIter = size(FFP_z, 2) - 2;
horizontalSignal = zeros(1, numIters*numSamplesPerIter+2);
FFP_z_fixed = zeros(1, numIters*numSamplesPerIter+2);
FFP_drive_fixed = zeros(1, numIters*numSamplesPerIter+2);
FFP_speed_fixed = zeros(1, numIters*numSamplesPerIter+2);
for k=1:numIters
    if k ~= numIters
        horizontalSignal(1, numSamplesPerIter*(k-1)+1:numSamplesPerIter*k) = horizontalSignal_mpi_mat(k, 1:end-2);
        FFP_z_fixed(1, numSamplesPerIter*(k-1)+1:numSamplesPerIter*k) = FFP_z(k, 1:end-2);
        FFP_drive_fixed(1, numSamplesPerIter*(k-1)+1:numSamplesPerIter*k) = FFP_drive(k, 1:end-2);
    else 
        horizontalSignal(1, numSamplesPerIter*(k-1)+1:numSamplesPerIter*k+2) = horizontalSignal_mpi_mat(k, :);
        FFP_z_fixed(1, numSamplesPerIter*(k-1)+1:numSamplesPerIter*k+2) = FFP_z(k, :);
        FFP_drive_fixed(1, numSamplesPerIter*(k-1)+1:numSamplesPerIter*k) = FFP_drive(k, 1:end-2);
    end
end

sigmoidParams = [3, 0];
t = (0:length(FFP_drive_fixed)-1)/MPIparams.fs;
drive_base = cos(2*pi*MPIparams.f_drive*t); % drive field movement
drive = sigmf(drive_base,sigmoidParams)-0.5;
drive = MPIparams.driveMag*drive/max(drive);

drive_der = -2*pi*MPIparams.f_drive*sin(2*pi*MPIparams.f_drive*t); % drive field derivative
drive_rep = MPIparams.driveMag*drive_der.*sigmf(drive_base,sigmoidParams).*(1-sigmf(drive_base,sigmoidParams));

L = length(drive_rep);
f = (0:L-1)*(MPIparams.fs)/L;
f_step = MPIparams.fs/L;
Bw_idx = round(20/f_step/2);

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
[filteredSignal, MPIparams] = extractData(horizontalSignal, MPIparams, 2e6, 20, 11);
[drive_rep_signal, MPIparams] = extractData(drive_rep, MPIparams, 2e6, 20, 11);

sig_fft = abs(fft(filteredSignal));
int_fft = abs(fft(drive_rep_signal));
sigPow = max(sig_fft(idx));
interferencePow = max(int_fft(idx));



sir = inf; snr = inf;
filteredSignal = filteredSignal;
interference = drive_rep_signal*sigPow/interferencePow/sir;

[signal, ~] = awgnInterference(filteredSignal, MPIparams, snr, inf);

plot(f, abs(fft(signal))); hold on; plot(f, abs(fft(interference)))

L = 800e3-648e3;
f = (0:L-1)*(MPIparams.fs)/L;
figure;
plot(f, abs(fft(signal(648e3+1:800e3)))); hold on; plot(f, abs(fft(interference(648e3+1:800e3))))


% start processing
start_idx = 1; % start_idx+val_idx-3; % corrected start idx "-2" is magical dont ask
pfov = (length(horizontalSignal)-2)/200; % total number of pfovs
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
for pfovIdx=900:915
    tic
    % sliding window
    % sig_idx = (1:(estimationParams.numSamplesPerIter-2)+2)+(estimationParams.numSample-2)*(pfovIdx-1); 
    % non-sliding window    
    sig_idx = (1:(estimationParams.numSamplesPerIter-2)+2)+(estimationParams.numSamplesPerIter-2)*(pfovIdx-1);
    FFP_z_part = FFP_z_fixed(sig_idx) - FFP_drive_fixed(sig_idx);
    robot_arm_pos(count) = FFP_z_part(floor(end/2)); 
    
    [tau_est_frequency(count), tau_est_linear(count), tau_lin_weighted(count)] = estimateTau_func(MPIparams, estimationParams, partial_signal_interp(sig_idx), 5);
%     err(count) = calculateError(MPIparams, estimationParams, partial_signal_interp(sig_idx), tau_est_linear(count));
    toc
    count = count + 1;
end

tau_estimation = zeros(1, length(tau_est_linear));
tau_estimation((robot_arm_pos <= SPIOparams.dz*200+SPIOparams.dz*10) & (robot_arm_pos >= SPIOparams.dz*200-SPIOparams.dz*10)) = ...
  tau_est_linear((robot_arm_pos <= SPIOparams.dz*200+SPIOparams.dz*10) & (robot_arm_pos >= SPIOparams.dz*200-SPIOparams.dz*10));
tau_estimation((robot_arm_pos <= -SPIOparams.dz*200+SPIOparams.dz*10) & (robot_arm_pos >= -SPIOparams.dz*200-SPIOparams.dz*10)) = ...
    tau_est_linear((robot_arm_pos <= -SPIOparams.dz*200+SPIOparams.dz*10) & (robot_arm_pos >= -SPIOparams.dz*200-SPIOparams.dz*10));
tau_estimation((robot_arm_pos <= SPIOparams.dz*10) & (robot_arm_pos >= -SPIOparams.dz*10)) = ...
    tau_est_linear((robot_arm_pos <= SPIOparams.dz*10) & (robot_arm_pos >= -SPIOparams.dz*10));

% tau_estimation((robot_arm_pos <= SPIOparams.dz*70+SPIOparams.dz*10) & (robot_arm_pos >= SPIOparams.dz*70-SPIOparams.dz*10)) = ...
%     tau_lin_weighted((robot_arm_pos <= SPIOparams.dz*70+SPIOparams.dz*10) & (robot_arm_pos >= SPIOparams.dz*70-SPIOparams.dz*10));
% tau_estimation((robot_arm_pos <= -SPIOparams.dz*70+SPIOparams.dz*10) & (robot_arm_pos >= -SPIOparams.dz*70-SPIOparams.dz*10)) = ...
%     tau_lin_weighted((robot_arm_pos <= -SPIOparams.dz*70+SPIOparams.dz*10) & (robot_arm_pos >= -SPIOparams.dz*70-SPIOparams.dz*10));
 
close all
figure; plot(robot_arm_pos*100, tau_estimation*1e6); axis tight, xlabel('z (cm)'); ylim([0, 5]);
ylabel('\tau (\mu s)', 'fontsize', 14); 
% title({'Cropped', ['Peak SIR = ' num2str(sir), ', Peak SNR = ' num2str(snr)]})
title({'Cropped Weighted Linear Estimation, 10 mm between Particles, ', ...
    ['3^{rd}-Harmonic SIR = ', num2str(sir), ', SNR = ' num2str(snr), ', R_s = ' num2str(MPIparams.slewRate), ' T/s, SIG = [10, 0]']})
figure; hold on;
% plot(robot_arm_pos*100, tau_est_frequency*1e6); 
% plot(robot_arm_pos*100, tau_est_linear*1e6);
plot(robot_arm_pos*100, tau_lin_weighted*1e6); axis tight;
xlabel('z (cm)'); ylabel('\tau (\mu s)', 'fontsize', 14); ylim([0, 10]); 
title({'Comparison Estimation, 10 mm between Particles, ', ...
    ['3^{rd}-Harmonic SIR = ', num2str(sir), ', SNR = ' num2str(snr), ', R_s = ' num2str(MPIparams.slewRate), ' T/s, SIG = [10, 0]']})
% legend('Frequency', 'Linear', 'Linear Weighted')
% title(['Peak SIR = ' num2str(sir), ', Peak SNR = ' num2str(snr)])





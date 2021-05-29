function [signal_filtered, MPIparams] = extractData(signal, MPIparams, interp_fs, Bw, numHarmonics)
    interp_coeff = interp_fs/MPIparams.fs;
    t = (0:length(signal)-1)/MPIparams.fs; t_interp = (0:(length(signal)-2)*interp_coeff+1)/interp_fs;
    signal = interp1(t, signal, t_interp, 'spline'); % interpolate the signal to an integer multiple of the drive frequency
     
    
    % signal = signal(round(t(end)*9.7e6*0.125):round(t(end)*9.7e6*0.875));

    MPIparams.fs = interp_fs;
    L = length(signal);
    f = (0:L-1)*(MPIparams.fs)/L;
    f_step = MPIparams.fs/L;


    % remove first harmonic and above 15th harmonic filter
    [~, idx_1] = min(abs(f-MPIparams.f_drive));
    [~, idx_2] = min(abs(f-(MPIparams.fs-MPIparams.f_drive)));
    Bw_idx = round(Bw/f_step/2);
    filter = ones(1, L);
    filter([1:idx_1+Bw_idx, idx_2-Bw_idx:L]) = 0;
    signal_no_1_f = fft(signal).*filter;
    [~, idx_1] = min(abs(f-MPIparams.f_drive*15));
    [~, idx_2] = min(abs(f-(MPIparams.fs-MPIparams.f_drive*15)));
    filter = ones(1, L);
    filter(idx_1:idx_2) = 0;
    signal_no_1_f = signal_no_1_f.*filter;
    signal_no_1 = ifft(signal_no_1_f);

    L = length(signal_no_1);
    f = (0:L-1)*(MPIparams.fs)/L;
    f_step = MPIparams.fs/L;
    Bw_idx = round(Bw/f_step/2);

    % extract odd harmonics
    idx = [];
    for k=3:2:numHarmonics
        [~, idx_1] = min(abs(f-MPIparams.f_drive*k));
        [~, idx_2] = min(abs(f-(MPIparams.fs-MPIparams.f_drive*k)));
        idx_1 = idx_1-Bw_idx:idx_1+Bw_idx;
        idx_2 = idx_2-Bw_idx:idx_2+Bw_idx;
        idx = [idx idx_1 idx_2];
    end
    filter = zeros(1, L);
    filter(idx) = 1;
    signal_f = fft(signal).*filter;
    signal_filtered = ifft(signal_f); % return back to time domain


%     % extract odd harmonics
%     idx = [];
%     for k=2:4
%         [~, idx_1] = min(abs(f-MPIparams.f_drive*k));
%         [~, idx_2] = min(abs(f-(MPIparams.fs-MPIparams.f_drive*k)));
%         idx_1 = idx_1-Bw_idx:idx_1+Bw_idx;
%         idx_2 = idx_2-Bw_idx:idx_2+Bw_idx;
%         idx = [idx idx_1 idx_2];
%     end
%     filter = zeros(1, L);
%     filter(idx) = 1;
%     signal_f = fft(signal).*filter;
%     signal_image = ifft(signal_f); % return back to time domain

end
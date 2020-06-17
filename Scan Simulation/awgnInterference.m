function [signal, interference] = awgnInterference(signal, MPIparams, snr, sir)

    interference = [];
    signal_f = fft(signal);
    L = length(signal);
    f = (0:L-1)*(MPIparams.fs)/L;
    f_step = MPIparams.fs/L;
    Bw = 20; % harmonic bandwidth (Hz)
    Bw_idx = round(Bw/f_step/2);
    
    % add interference
    if ~isinf(sir)
        % extract odd harmonics and add interference
        interference = zeros(1, L);
        idx = [];
        for k=3:2:7
            [~, idx_1_c] = min(abs(f-MPIparams.f_drive*k));
            [~, idx_2_c] = min(abs(f-(MPIparams.fs-MPIparams.f_drive*k)));
            idx_1 = idx_1_c-Bw_idx:idx_1_c+Bw_idx;
            idx_2 = idx_2_c-Bw_idx:idx_2_c+Bw_idx;
            idx = [idx_1 idx_2];
            harmonicPower = max(abs(signal_f(idx)));
            interference([idx_1_c idx_2_c]) = harmonicPower/sir*exp(1i*2*pi*rand(1, 2));
        end
        signal = ifft(signal_f + interference); % return back to time domain
    end

    % add noise
    if ~isinf(snr)
        [sigPow, ~] = max(abs(signal));
        sigma = sigPow/snr;
        signal=signal+sigma*randn(1,length(signal));
    end
end
function signal_filtered = signalFrequencySelection(signal,fs,fd,num_harmonics, harmonics, bandwidth, signal_type)

    % if signal is already in frequency domain, do nothing if not take FFT
    switch signal_type
        case 'Time'
            signal = fft(signal, size(signal, 1));
        case 'Frequency'
            
    end
    N = (size(signal, 1)-2)/fs*fd;

    switch harmonics
        case 'all'
            I = [2:1:num_harmonics]*N;
            I = [round(I+1) round(length(signal)-fliplr(I)+1)];
        case 'odd'
            I = [3:2:num_harmonics]*N;
            I = [round(I+1) round(length(signal)-fliplr(I)+1)];
        case 'first'
            I = 1*N;
            I = [round(I+1) round(length(signal)-fliplr(I)+1)];
        case 'linear'
            I = [1:1:num_harmonics]*N;
            I = [round(I+1) round(length(signal)-fliplr(I)+1)];

            I = round(bsxfun(@plus, I', -bandwidth/2:bandwidth/2));
            I = transpose(I(:));
    %             I_orig = I;
    %             for k=1:length(I_orig)
    %                 I(I == I_orig(k)) = [];
    %             end
        case 'no_filter'
            I = 1:length(signal);
    end

    % removes noise from irrelevant frequencies
    fft_signal_filt = zeros(size(signal));
    fft_signal_filt(I, :) = signal(I, :);
    signal_filtered = ifft(fft_signal_filt, size(signal, 1), 'symmetric');
    
end
function [err] = calculateError(MPIparams, estimationParams, signal, tau_est)

    if tau_est < 0
        tau_est = 1e-10;
    end
    t = (0:1/MPIparams.fs:tau_est*15);
    r_t = 1/tau_est*exp(-t./tau_est);
    r_t = r_t/sum(r_t);
    r_f = fft(r_t, estimationParams.numSamplesPerIter);
    sig_f = fft(signal);
    sig = ifft(sig_f./r_f);
    
    pos = sig(2:end-MPIparams.fs/MPIparams.f_drive/2-2); 
    neg = sig(MPIparams.fs/MPIparams.f_drive/2+3:end-1);
    L = length(neg);
    f = (0:L-1)*(MPIparams.fs)/L-(MPIparams.fs)/2;
    f = fftshift(f);
    neg = real(ifft(fft(neg).*exp(1i*2*pi*estimationParams.del_t.*f)))*estimationParams.amp_t;
    err = (norm(pos+flip(neg), 2))^0.5;
       
end
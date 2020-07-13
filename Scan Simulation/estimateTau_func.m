function [tau_est_frequency, tau_est_linear, tau_lin_weighted] = estimateTau_func(MPIparams, estimationParams, signal, freqContribution)


    idx = 1:estimationParams.numSamplesPerIter;
    t_sig = idx/MPIparams.fs;

    idx_interp = (1:estimationParams.numSampleInterpolated);
    t_interp = idx_interp/MPIparams.fs/estimationParams.interp_coeff;

    sig = interp1(t_sig, signal,t_interp, 'linear', 'extrap');


    
    pos = sig((1:end-MPIparams.fs/MPIparams.f_drive/2*estimationParams.interp_coeff-2*estimationParams.interp_coeff)); 
    neg = sig(MPIparams.fs/MPIparams.f_drive/2*estimationParams.interp_coeff+(2*estimationParams.interp_coeff+1):end); 

    L = length(neg);
    f = (0:L-1)*(MPIparams.fs*estimationParams.interp_coeff)/L-(MPIparams.fs*estimationParams.interp_coeff)/2;
    f = fftshift(f);

    S1 = fft(pos);
    S2 = fft(neg).*exp(1i*2*pi*estimationParams.del_t.*f)*estimationParams.amp_t;

    sum_val = (S1+conj(S2));
    sub_val = (conj(S2)-S1);
    
    

    a = 1i*2*pi*f.*sub_val;
    b = sum_val;
    


    f_axis = (0:L-1)*(MPIparams.fs*estimationParams.interp_coeff)/L;

    sig_contribution = [];
    for cont_idx=3:2:freqContribution
        [~, idx_1] = min(abs(f_axis-MPIparams.f_drive*cont_idx));
        [~, idx_2] = min(abs(f_axis-(MPIparams.fs-MPIparams.f_drive*cont_idx)));
        sig_contribution = [sig_contribution idx_1+1];
    end


    a = transpose(real(a(sig_contribution)));
    b = transpose(real(b(sig_contribution)));

    tau_est_frequency = mean(real(b./a));
    tau_est_linear = (a'*a)^-1*a'*b; 

    W = diag(abs(S1(sig_contribution)).^2);
    tau_lin_weighted = (a'*W*a)^-1*(a'*W)*b;

end
function [del_t, a_t] = findDistortion(MPIparams)
    Bp = MPIparams.Bp;
    Gzz = MPIparams.Gzz;
    R = MPIparams.Rs(3)/MPIparams.Gzz;
    f_d = MPIparams.f_drive;

    arg = 2*pi*f_d;
    % taylor series expension of sine
    N = 5;
    p = zeros(1, N+1); 
    for k=0:N
        p(k+1) = ((-1)^k)/factorial(2*k+1)*(arg^(2*k+1)); 
    end
    p = p*MPIparams.driveMag;
    p = flip(p);
    p = upsample(p, 2); % turn into polynomial
    p(end-1) = p(end-1) + R; % add the extra R*del_t 
    p(end) = p(end) + R/(2*f_d);
    rp = roots(p); % find roots
    real_idx = imag(rp) == 0; % extract real roots
    real_roots = rp(real_idx);
    del_t = real_roots(end); % select the last one (last one is always the root we are looking for)

    a_t = (Bp/Gzz*2*pi*f_d*cos(2*pi*f_d*del_t) - R)/(Bp/Gzz*2*pi*f_d + R); % amplitude distortion

end
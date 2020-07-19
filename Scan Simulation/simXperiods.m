function xPeriodsSim = simXperiods(MPIparams)

    xSpeed = 1/(MPIparams.Gxx*(1/MPIparams.Rs(1)));
    xPeriods = MPIparams.numTrianglePeriods*2*MPIparams.FOV_x*MPIparams.f_drive/xSpeed;
    t = (0:xPeriods-1)/MPIparams.f_drive;
    p = 1/(MPIparams.numTrianglePeriods/MPIparams.time);
    xPeriodsPos = MPIparams.FOV_x*(2*abs(2*(t/p - floor(t/p + 0.5)))-1)/2; % robot arm movement in x direction w.r.t. time

    xPeriodsSim = [];
    for k=1:length(MPIparams.traversedFOVx)/2
        val = [MPIparams.traversedFOVx((k-1)*2+1) MPIparams.traversedFOVx(k*2)];
        idx = ((xPeriodsPos >= val(1)) & (xPeriodsPos <= val(2)));
        xPeriodsSim = [xPeriodsSim, find(idx==1)];
    end

end
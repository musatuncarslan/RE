function xPeriodsSim = simXperiods(MPIparams)


    xPeriods = MPIparams.zPeriods/MPIparams.numTrianglePeriods/2; % number of periods per line segment in x-axis
    xPeriodsSim = [];
    for k=1:MPIparams.numTrianglePeriods*2
        t = ((k-1)*xPeriods:k*xPeriods-1)/MPIparams.f_drive;
        p = 1/(MPIparams.numTrianglePeriods/MPIparams.time);
        xPeriodsPos = MPIparams.FOV_x*(2*abs(2*(t/p - floor(t/p + 0.5)))-1)/2; % robot arm movement in x direction w.r.t. time

        [~, closestPoint] = min(abs(xPeriodsPos-MPIparams.traversedFOVx'), [], 2);
        if mod(k, 2); else; closestPoint = flip(closestPoint); end
        closestPoint = closestPoint + xPeriods*(k-1);
        for l=1:length(closestPoint)/2
            xPeriodsSim = [xPeriodsSim, closestPoint((l-1)*2+1):closestPoint(l*2)];
        end
    end
    
%     Simparams.xPeriodsSim = xPeriodsSim;

end
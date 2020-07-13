function zPeriodsSim = simZperiods(MPIparams)

    zSpeed = 1/(MPIparams.Gzz*(1/MPIparams.Rsz));
    zPeriodsPos = (0:MPIparams.zPeriods-1)/MPIparams.f_drive*zSpeed - MPIparams.FOV_z/2;
    [~, closestPoint] = min(abs(zPeriodsPos-MPIparams.traversedFOVz'), [], 2);
    zPeriodsSim = [];
    for k=1:length(closestPoint)/2
        zPeriodsSim = [zPeriodsSim, closestPoint((k-1)*2+1):closestPoint(k*2)];
    end
%     Simparams.zPeriodsSim = zPeriodsSim;

end

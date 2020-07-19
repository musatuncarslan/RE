function zPeriodsSim = simZperiods(MPIparams)

    zSpeed = 1/(MPIparams.Gzz*(1/MPIparams.Rs(3)));
    zPeriodsPos = (0:MPIparams.zPeriods-1)/MPIparams.f_drive*zSpeed - MPIparams.FOV_z/2;
    zPeriodsSim = [];
    for k=1:length(MPIparams.traversedFOVz)/2
        val = [MPIparams.traversedFOVz((k-1)*2+1) MPIparams.traversedFOVz(k*2)];
        idx = ((zPeriodsPos >= val(1)) & (zPeriodsPos <= val(2)));
        zPeriodsSim = [zPeriodsSim, find(idx==1)];
    end
end

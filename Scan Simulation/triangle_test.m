clear all
close all
clc

format long;

gpudev = gpuDevice(1); % get the GPU device

% parameters
Physicsparams = setPhysicsParams(); % physics parameters
MPIparams = setMPIParams(Physicsparams, 'complex_rastered', 0.1); % MPI machine parameters
SPIOparams = setSPIOParams(Physicsparams, 512, 2e-6); % SPIO parameters
Simparams = setSimulationParams(MPIparams, Physicsparams); % Simulation parameters

particleNo = 1;

% preprocessing in order to find all the angles to calculate the necessary
% PSFs
uniqueAngle_sim = [];
count = 1;
for k=1   
    t = gpuArray(((Simparams.simPeriods(1)-1)*Simparams.samplePerPeriod+1:(Simparams.simPeriods(end)*Simparams.samplePerPeriod+2*Simparams.downsample+1))/Simparams.fs_phsy);
    FFPparams = generateFFP(gpudev, t, MPIparams, Simparams, [3, 0]); 
    uniqueAngle_sim = [uniqueAngle_sim, FFPparams.FFP_uniqueAngle];
    count = count + 1;
end
angleVec = unique(uniqueAngle_sim);

x_axis = (-SPIOparams.image_FOV_x/2:SPIOparams.dx:SPIOparams.image_FOV_x/2-SPIOparams.dx);
z_axis = (-SPIOparams.image_FOV_z/2:SPIOparams.dz:SPIOparams.image_FOV_z/2-SPIOparams.dz);

for k=1:length(SPIOparams.diameter)
    surf(x_axis, z_axis, SPIOparams.SPIOdistribution(:,:,k)); shading interp
    hold on;
end
plot(FFPparams.FFP_x, FFPparams.FFP_z)
view(2)
xlabel('x-axis'); ylabel('z-axis')
xlim([-MPIparams.FOV_x/2 MPIparams.FOV_x/2]); ylim([-MPIparams.FOV_z/2 MPIparams.FOV_z/2])

wait(gpudev); clear FFPparams uniqueAngle_sim t;

div = divisors(length(Simparams.simPeriods));
big = div(end);
if isempty(div(div > big)) == 1
    div = div(end);
else
    div = div(div > big); 
    div = div(1);
end
divL = div;
div = length(Simparams.simPeriods)/div;

signal = zeros(div, Simparams.samplePerPeriod*length(Simparams.simPeriods)/Simparams.downsample/div+2);
FFP_x = zeros(div, Simparams.samplePerPeriod*length(Simparams.simPeriods)/Simparams.downsample/div+2);
FFP_z = zeros(div, Simparams.samplePerPeriod*length(Simparams.simPeriods)/Simparams.downsample/div+2);
drive = zeros(div, Simparams.samplePerPeriod*length(Simparams.simPeriods)/Simparams.downsample/div+2);
for particleNo = 1:length(SPIOparams.diameter)
    % efficient functions, uses HDD instead of RAM, should not give any
    % memory errors, uses HDF5 file structure
    tic
    generatePSFe(gpudev, MPIparams, SPIOparams, angleVec, particleNo, 50); 
    toc

    % generate the MPI signal
    for k=1:div
        tic
        idx = [(k-1)*divL+1, k*divL];
        simIdx = ((Simparams.simPeriods(idx(1))-1)*Simparams.samplePerPeriod+1:(Simparams.simPeriods(idx(end))*Simparams.samplePerPeriod+2*Simparams.downsample+1));
        t = gpuArray(simIdx/Simparams.fs_phsy);
        FFPparams = generateFFP(gpudev, t, MPIparams, Simparams, [3, 0]);   
        wait(gpudev); clear t;
        [signals_sep, SPIOparams] = generateSe(gpudev, FFPparams, MPIparams, SPIOparams, Simparams, particleNo, 50);
        toc
        signal(k, :) = signal(k,:) + gather(signals_sep.horizontalSignal); 
        wait(gpudev); clear FFPparams t signals_sep;
        
        if particleNo == 1
            FFP_x(k, :) = gather(FFP_params.FFP_x(1:Simparams.downsample:end-1));
            FFP_z(k, :) = gather(FFP_params.FFP_z(1:Simparams.downsample:end-1));
            drive(k, :) = gather(FFP_params.FFP_drive(1:Simparams.downsample:end-1));
        end
    end


end

% delete signal file inside the temporary folder if it exists
tempList = dir('temp');
if length(tempList) > 3
    delete('temp\signal.h5')
end

sigL = size(signal);
h5create(['./temp/','signal.h5'],'/signal',sigL);
h5create(['./temp/','signal.h5'],'/FFP_x',sigL);
h5create(['./temp/','signal.h5'],'/FFP_z',sigL);
h5create(['./temp/','signal.h5'],'/drive',sigL);

h5write('./temp/signal.h5', '/signal', signal, [1 1], sigL);
h5write('./temp/signal.h5', '/FFP_x', FFP_x, [1 1], sigL);
h5write('./temp/signal.h5', '/FFP_z', FFP_z, [1 1], sigL);
h5write('./temp/signal.h5', '/drive', drive, [1 1], sigL);
% tic
% % regular PSF generation, uses RAM instead of HDD, may give memory errors
% % if angleVec is too long (i.e. angle resolution is high)
% [colinearPSF1,transversePSF1] = generatePSF(MPIparams, SPIOparams, FFPparams);
% [colinearIMG, transverseIMG] = generatePSFImages(SPIOparams, colinearPSF1, transversePSF1);
% toc

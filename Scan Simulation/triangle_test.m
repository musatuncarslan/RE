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

uniqueAngle_sim = [];
count = 1;
for k=Simparams.simPeriods   
    t = (((k-1)*Simparams.samplePerPeriod+1:(k*Simparams.samplePerPeriod+2*Simparams.downsample+1))/Simparams.fs_phsy);
    FFPparams = generateFFP(t, MPIparams, Simparams, [3, 0]);  
    uniqueAngle_sim = [uniqueAngle_sim, FFPparams.FFP_uniqueAngle];
    count = count + 1;
end
angleVec = unique(uniqueAngle_sim);

particleNo = 1;

% tic
% % efficient functions, uses HDD instead of RAM, should not give any
% % memory errors, uses HDF5 file structure
% generatePSFe(gpudev, MPIparams, SPIOparams, angleVec, particleNo, 50); 
% generatePSFIe(gpudev, SPIOparams, particleNo, 50)
% toc


for k=Simparams.simPeriods
    tic
    t = gpuArray(((k-1)*Simparams.samplePerPeriod+1:(k*Simparams.samplePerPeriod+2*Simparams.downsample+1))/Simparams.fs_phsy);
    FFPparams = generateFFP(t, MPIparams, Simparams, [3, 0]);   
    [signals_sep, SPIOparams] = generateSe(angleVec, FFPparams, MPIparams, SPIOparams, Simparams, particleNo, 10);
    toc
end


% tic
% % regular PSF generation, uses RAM instead of HDD, may give memory errors
% % if angleVec is too long (i.e. angle resolution is high)
% [colinearPSF1,transversePSF1] = generatePSF(MPIparams, SPIOparams, FFPparams);
% [colinearIMG, transverseIMG] = generatePSFImages(SPIOparams, colinearPSF1, transversePSF1);
% toc

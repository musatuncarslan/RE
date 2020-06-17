clear all
close all
clc

format long;

addpath('./parameterFiles')
addpath('./signalGenerationFiles')

% parameters
Physicsparams = setPhysicsParams(); % physics parameters
MPIparams = setMPIParams(Physicsparams, 0.1); % MPI machine parameters

SPIOparams = setSPIOParams(Physicsparams, 512, 2e-6); % SPIO parameters
Simparams = setSimulationParams(MPIparams, Physicsparams); % Simulation parameters

filePath = ''; % 'C:\Users\Orion\Google Drive\Phd\signalSaveFolder\';
fileName = ['signal_' datestr(now,'dd-mmmm-yyyy_HH-MM-SS') '.mat'];
fullFileName = [filePath fileName];
[fileObj, signal_size, ~] = setFileParams(fullFileName, Physicsparams, MPIparams, SPIOparams, Simparams);

chunk = signal_size;
count = 1;
sig = [];
for k=Simparams.startIter:Simparams.endIter    
    tic
    t = gpuArray(((k-1)*Simparams.numSamplesPerIter+1:(k*Simparams.numSamplesPerIter+2*Simparams.downsample+1))/Physicsparams.fs);
    FFPparams = generateFFP(t, MPIparams, Simparams, [3, 0]);  
    % in a single period, necessary angles will probably be gathered so
    % generating PSFs and corresponding images for a single period should
    % be enough
    if count == 1
        % calculate colinear and trasnverse PSF(s) for each unique angle
        [colinearPSF, transversePSF, X, Z] = generatePSF(MPIparams, SPIOparams, FFPparams);

        % pad the empty parts of the image so that it has same size with PSF(s),
        % then compute colinear and transverse images for each unique angle
        [colinearIMG, transverseIMG] = generatePSFImages(SPIOparams, colinearPSF, transversePSF);
    end
    % generate MPI signals
    [signals_sep, SPIOparams] = generateSignals(colinearIMG, transverseIMG, FFPparams, MPIparams, SPIOparams, Simparams);
      
    
    signals = struct;
    signals.horizontalSignal = zeros(1, size(signals_sep, 2));
    signals.verticalSignal = zeros(1, size(signals_sep, 2));
    for l=1:length(SPIOparams.diameter)
        signals.horizontalSignal = signals.horizontalSignal + real(signals_sep.horizontalSignal(l, :));
        signals.verticalSignal = signals.verticalSignal + real(signals_sep.verticalSignal(l, :));
    end
    % save signals to the file
    fileObj.horizontalSignal_mpi_mat(count,1:signal_size) = signals.horizontalSignal;
    fileObj.verticalSignal_mpi_mat(count,1:signal_size) = signals.verticalSignal;
    fileObj.FFP_x(count,1:signal_size) = FFPparams.FFP_xToFile;
    fileObj.FFP_z(count,1:signal_size) = FFPparams.FFP_zToFile;
    fileObj.FFP_drive(count,1:signal_size) = FFPparams.FFP_driveToFile;
    fileObj.FFP_speed(count,1:signal_size) = FFPparams.FFP_speedToFile;
    fileObj.FFP_angle(count,1:signal_size) = FFPparams.FFP_angleToFile;
    
    endTime(count) = toc;
    fprintf('%2.2f %% done in %2.4f sec.\n', 100*chunk/signal_size/Simparams.numIters, endTime(count));
    chunk = chunk + signal_size;
    count = count + 1;
    
    
    
end

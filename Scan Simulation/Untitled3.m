clear all
close all
clc

format long;
format compact;

gpudev = gpuDevice(1); % get the GPU device

rs_z = 0.1; % [0.1, 0.5, round(linspace(1, 10, 10), 2)];
rs_x = [0.1, 0.5, round(linspace(1, 10, 10), 2)];

for rs_z_idx = rs_z
   for rs_x_idx = rs_x
        % parameters
        Physicsparams = setPhysicsParams(); % physics parameters
        MPIparams = setMPIParams(Physicsparams, [rs_x_idx, 0, rs_z_idx]); % MPI machine parameters
        SPIOparams = setSPIOParams(Physicsparams, 512, 2e-6); % SPIO parameters
        [Simparams, MPIparams] = setSimulationParams(MPIparams, Physicsparams); % Simulation parameters

        numPeriods = 4; % seconds

        uniqueAngle_sim = [];
        t = gpuArray((1:Simparams.samplePerPeriod*numPeriods+2*Simparams.downsample+1)/Simparams.fs_phsy);
        FFPparams = generateFFP(gpudev, t, MPIparams, 'test'); 
        uniqueAngle_sim = [uniqueAngle_sim, FFPparams.FFP_uniqueAngle];
        FFP_x = gather(FFPparams.FFP_x(1:Simparams.downsample:end-1));
        FFP_z = gather(FFPparams.FFP_z(1:Simparams.downsample:end-1));
        angleVec = unique(uniqueAngle_sim);

        figure; hold on;
        x_axis = (-SPIOparams.image_FOV_x/2:SPIOparams.dx:SPIOparams.image_FOV_x/2-SPIOparams.dx);
        z_axis = (-SPIOparams.image_FOV_z/2:SPIOparams.dz:SPIOparams.image_FOV_z/2-SPIOparams.dz);
        for k=1:length(SPIOparams.diameter)
            surf(x_axis*100, z_axis*100, SPIOparams.SPIOdistribution(:,:,k)); shading interp
            hold on;
        end
        plot(FFP_x*100, FFP_z*100, 'color', [0 0.4470 0.7410])
        view(2)
        xlabel('x-axis (cm)'); ylabel('z-axis (cm)')
        % xlim([-MPIparams.FOV_x/2 MPIparams.FOV_x/2]); ylim([-MPIparams.FOV_z/2 MPIparams.FOV_z/2])

        wait(gpudev); clear FFPparams uniqueAngle_sim t;

        % actual signal generation
        signal = zeros(Simparams.divNum, numPeriods*Simparams.samplePerPeriod/Simparams.downsample+2);
        sigL = size(signal);
        FFP_x = zeros(sigL);
        FFP_z = zeros(sigL);
        drive = zeros(sigL);
        div = Simparams.div;
        for particleNo = 1:length(SPIOparams.diameter)
            % efficient functions, uses HDD instead of RAM, should not give any
            % memory errors, uses HDF5 file structure
            generatePSFe(gpudev, MPIparams, SPIOparams, angleVec, particleNo, 40); 

            % generate the MPI signal
            for k=1:1
                t = gpuArray((1:Simparams.samplePerPeriod*numPeriods+2*Simparams.downsample+1)/Simparams.fs_phsy);
                FFPparams = generateFFP(gpudev, t, MPIparams, 'test'); 
                wait(gpudev); clear t;
                [signals_sep, SPIOparams] = generateSe(gpudev, FFPparams, MPIparams, SPIOparams, Simparams, particleNo, 50); % simulate signal

                signal(k, 1:length(signals_sep.horizontalSignal)) = signal(k,1:length(signals_sep.horizontalSignal)) + gather(signals_sep.horizontalSignal); 
                wait(gpudev); clear t;

                if particleNo == 1
                    FFP_x(k, 1:length(signals_sep.horizontalSignal)) = gather(FFPparams.FFP_x(1:Simparams.downsample:end-1));
                    FFP_z(k, 1:length(signals_sep.horizontalSignal)) = gather(FFPparams.FFP_z(1:Simparams.downsample:end-1));
                    drive(k, 1:length(signals_sep.horizontalSignal)) = gather(FFPparams.drive(1:Simparams.downsample:end-1));
                end
                wait(gpudev); clear signals_sep FFPparams;
            end
        end

        % delete signal file inside the temporary folder if it exists
        tempList = dir('temp');
        if length(tempList) > 3
            delete('temp\signal.h5')
        end


        h5create(['./temp/','signal.h5'],'/signal',sigL);
        h5create(['./temp/','signal.h5'],'/FFP_x',sigL);
        h5create(['./temp/','signal.h5'],'/FFP_z',sigL);
        h5create(['./temp/','signal.h5'],'/drive',sigL);

        h5write('./temp/signal.h5', '/signal', signal, [1 1], sigL);
        h5write('./temp/signal.h5', '/FFP_x', FFP_x, [1 1], sigL);
        h5write('./temp/signal.h5', '/FFP_z', FFP_z, [1 1], sigL);
        h5write('./temp/signal.h5', '/drive', drive, [1 1], sigL);

        clearvars -except Physicsparams MPIparams SPIOparams Simparams rs_z_idx rs_z rs_x_idx rs_x gpudev
        estimateTau_test
        close all
   end
end


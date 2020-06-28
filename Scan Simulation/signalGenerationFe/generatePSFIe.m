function generatePSFIe(gpudev, SPIOparams, particleNo, angleVec, div)

    info = h5info('./temp/PSF.h5');
    dataSize = info.Datasets(1).Dataspace.Size;
    zL = dataSize(1); xL = dataSize(2);
    % divide data into chunks
    numAngle = dataSize(3);

    numIters = floor(numAngle/div);
    if ((numAngle/div) > numIters)
        numIters = numIters + 1;
    end
    iterVec = [div*ones(1, numIters-1), mod(numAngle, div)];
    idxVec = [0 cumsum(iterVec)];

    tempDistribution = SPIOparams.SPIOdistribution(:,:,particleNo);

    % pad the distribution with 0s so it has same size with the PSFs
    img_size = size(tempDistribution); psf_size = [dataSize(1), dataSize(2)];
    if img_size(1) < psf_size(1)
       tempDistribution = padarray(tempDistribution,[floor((psf_size(1)-img_size(1))/2), 0],0,'both');
       img_size = size(tempDistribution);
       tempDistribution = padarray(tempDistribution,[psf_size(1)-img_size(1), 0],0,'post');
    end
    if img_size(2) < psf_size(2)
       tempDistribution = padarray(tempDistribution,[0, floor((psf_size(2)-img_size(2))/2)],0,'both');
       img_size = size(tempDistribution);
       tempDistribution = padarray(tempDistribution,[0, psf_size(2)-img_size(2)],0,'post');
    end
    % multiplication in frequency domain is faster than convolution in
    % time domain
    imgF = gpuArray(fft2(tempDistribution, dataSize(1), dataSize(2)));

    % create new datasets for colinear and transverse images
    h5create(['./temp/','IMG.h5'],'/colinearIMG',[zL xL numAngle]);
    h5create(['./temp/','IMG.h5'],'/transverseIMG',[zL xL numAngle]);
    h5create(['./temp/','IMG.h5'],'/angleVec',[1 length(angleVec)]);
    for l=1:numIters
        colinearPSF_f = fft2(gpuArray(h5read('./temp/PSF.h5','/colinearPSF', [1 1 idxVec(l)+1], [zL xL iterVec(l)])));
        colIMG = real(ifft2(imgF.*colinearPSF_f)); wait(gpudev); clear colinearPSF_f;
        colIMG_s = fftshift(fftshift(colIMG, 1), 2); wait(gpudev); clear colIMG;
        h5write('./temp/IMG.h5', '/colinearIMG', gather(colIMG_s), [1 1 idxVec(l)+1], [1201 1001 iterVec(l)]); wait(gpudev); clear colIMG_s;

        transversePSF_f = fft2(gpuArray(h5read('./temp/PSF.h5','/transversePSF', [1 1 idxVec(l)+1], [zL xL iterVec(l)])));
        tranIMG = real(ifft2(imgF.*transversePSF_f)); wait(gpudev); clear transversePSF_f;
        tranIMG_s = fftshift(fftshift(tranIMG, 1), 2); wait(gpudev); clear tranIMG;
        h5write('./temp/IMG.h5', '/transverseIMG', gather(tranIMG_s), [1 1 idxVec(l)+1], [1201 1001 iterVec(l)]); wait(gpudev); clear tranIMG_s;
    end
    wait(gpudev); clear tranIMG_s;
    h5write('./temp/IMG.h5', '/angleVec', gather(angleVec), [1 1], [1 length(angleVec)]); % write angle vector to file for ease of access later on

    % delete the unnecessary PSF file to clear up space
    tempList = dir('temp');
    if length(tempList) > 2
        delete('temp\PSF.h5');
    end

end
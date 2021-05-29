function SPIOparams = setSPIOParams(Physicsparams, pos, tau_val)
    SPIOparams = struct;
    
    SPIOparams.diameter = [25]; % (nm)
    SPIOparams.tau = [4e-6, 0, 3e-6]; % (S)
    
    thc = 0;
    hei = 0;
    
    % spio distribution in 2D
    SPIOparams.SPIOdistribution = zeros(512, 512, length(SPIOparams.diameter)); 
    SPIOparams.SPIOdistribution((256-hei:256+hei), (256-thc:256+thc), 1) = 1;
%     for k=-5:1:5
%         SPIOparams.SPIOdistribution(251-k+50, round((251-thc:251+thc)+10*k), 1) = 1;
%     end
%     SPIOparams.SPIOdistribution((350-hei:350+hei), (257-thc:257+thc), 2) = 1.5;  
%     SPIOparams.SPIOdistribution = imrotate(SPIOparams.SPIOdistribution, 22);
    
%     SPIOparams.SPIOdistribution((251-20:251+20)+10, (301-20:301+20), 3) = 1;

    % extent of SPIO distribution 
    SPIOparams.image_FOV_x = 0.05; % (m) 
    SPIOparams.image_FOV_z = 0.05; % (m) 
    SPIOparams.dx = SPIOparams.image_FOV_x/size(SPIOparams.SPIOdistribution,1);  % distance between each pixel (m)
    SPIOparams.dz = SPIOparams.image_FOV_z/size(SPIOparams.SPIOdistribution,2);  

    SPIOparams.r_t = cell(1, length(SPIOparams.diameter));
    for k=1:length(SPIOparams.diameter)
        if (SPIOparams.tau(k) == 0)
            SPIOparams.r_t{k} = 1;
            SPIOparams.horizontalPrev{k} = 0;
            SPIOparams.verticalPrev{k} = 0;            
        else
            t = (0:1/Physicsparams.fs:SPIOparams.tau(k)*15);
            SPIOparams.r_t{k} = 1/SPIOparams.tau(k)*exp(-t./SPIOparams.tau(k));
            SPIOparams.horizontalPrev{k} = zeros(1, length(t));
            SPIOparams.verticalPrev{k} = zeros(1, length(t));       
        end
    end
end
clear all
close all
clc

fs = 0.5e6;
f_drive = 10.8e3;
FOV_x = 0.05;
FOV_z = 0.06;
time = 1.43; % time to traverse whole FOV
numTrianglePeriods = 3; % number of triangle zig-zags in whole movement

fs = fs-mod(fs, f_drive);
numPeriods = time*f_drive; % number of drive periods during whole FOV scan (this is hopefully an integer)
numPeriodsPerTriangle = numPeriods/numTrianglePeriods; % number of drive periods per triangle periods 
                                                       % (regarding to the movement in x-direction, 
                                                       % this is also hopefully an even integer, if necessary adjust time so that it is an even integer)

if (floor(numPeriods)~=numPeriods)
    error('Error. Make sure that time*f_drive is an integer. \n\n \t You can slightly adjust time for this.')
end
if (floor(numPeriodsPerTriangle)~=numPeriodsPerTriangle)
    error('Error. Make sure that time*f_drive/numTrianglePeriods is an integer. \n\n \t You can slightly adjust time for this.')
end
if mod(numPeriodsPerTriangle, 2) == 1
    error('Error. Make sure that time*f_drive/numTrianglePeriods is divisible by 2. \n\n \t You can slightly adjust time for this.')
end     
                                                       
samplePerPeriods = fs*time/numPeriods; % number of samples per period

robotSpeed_z = FOV_z/time; % robot arm movement speed in z direction (m/s)
robotSpeed_x = FOV_x*numTrianglePeriods*2/time; % robot arm movement speed in x direction (m/s)

traversedFOVz = [-0.03 0.03]; % z-axis FOV that is of interest (m)
traversedFOVx = [-0.025 -0.02]; % x-axis FOV that is of interest (m)

xi = robotSpeed_x*time/numPeriods;

x11 = round((traversedFOVx(1)+FOV_x/2)/xi);
x12 = round((traversedFOVx(2)+FOV_x/2)/xi)-1;

x21 = round((-traversedFOVx(2)+FOV_x/2)/xi)+numPeriodsPerTriangle/2;
x22 = round((-traversedFOVx(1)+FOV_x/2)/xi)+numPeriodsPerTriangle/2-1;

L = ((x12-x11)+(x22-x21))+2;
simPeriodsx = zeros(1, L*numTrianglePeriods);
for k=1:numTrianglePeriods
    simPeriodsx((k-1)*L+1:k*L) = ([x11:x12, x21:x22])+(k-1)*numPeriodsPerTriangle; % simulated periods
end

numTrianglePeriods = numTrianglePeriods/time;
p = 1/numTrianglePeriods;

x_partial = []; % zeros(length(simPeriods)*samplePerPeriods, 1);
z_partial = []; % zeros(length(simPeriods)*samplePerPeriods, 1);
for k=simPeriodsx
    t = (k*samplePerPeriods:((k+1)*samplePerPeriods-1))/fs; 
    x_partial = [x_partial FOV_x*(2*abs(2*(t/p - floor(t/p + 0.5)))-1)/2]; % robot arm movement in x direction w.r.t. time
    z_partial = [z_partial FOV_z/time*t - FOV_z/2]; % robot arm movement in z direction w.r.t. time
end


t = (0:fs*time-1)/fs;
x = FOV_x*(2*abs(2*(t/p - floor(t/p + 0.5)))-1)/2;
z = FOV_z/time*t - FOV_z/2;

figure; plot(z_partial, x_partial, '*'); hold on;
plot(z, x, '.-');
xlabel('z-axis'); ylabel('x-axis')
xlim([-0.03 0.03]); ylim([-0.025 0.025])

legend('Partial Trajectory', 'Complete Trajectory')
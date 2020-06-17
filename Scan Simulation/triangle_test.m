clear all
close all
clc

fs = 0.5e6;
f_drive = 50e3;
FOV_x = 0.05;
FOV_z = 0.06;
time = 2; % time to traverse whole FOV
numTrianglePeriods = 1; % number of triangle zig-zags in whole movement


numPeriods = time*f_drive; % number of drive periods during whole FOV scan
numPeriodsPerTriangle = numPeriods/numTrianglePeriods; % number of drive periods per triangle periods (regarding to the movement in x-direction)
samplePerPeriods = fs*time/numPeriods; % number of samples per period

robotSpeed_z = FOV_z/time; % robot arm movement speed in z direction (m/s)
robotSpeed_x = FOV_x*numTrianglePeriods*2/time; % robot arm movement speed in x direction (m/s)

traversedFOVz = [-0.03 0.03]; % z-axis FOV that is of interest (m)
traversedFOVx = [0.02 0.025]; % x-axis FOV that is of interest (m)

xi = robotSpeed_x*time/numPeriods;

k11 = round((traversedFOVx(1)+FOV_x/2)/xi);
k12 = round((traversedFOVx(2)+FOV_x/2)/xi)-1;

k21 = round((-traversedFOVx(2)+FOV_x/2)/xi)+numPeriodsPerTriangle/2;
k22 = round((-traversedFOVx(1)+FOV_x/2)/xi)+numPeriodsPerTriangle/2-1;


simPeriods = [k11:k12, k21:k22]; % simulated periods

numTrianglePeriods = numTrianglePeriods/2;
p = 1/numTrianglePeriods;

x_partial = []; % zeros(length(simPeriods)*samplePerPeriods, 1);
z_partial = []; % zeros(length(simPeriods)*samplePerPeriods, 1);
for k=simPeriods
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
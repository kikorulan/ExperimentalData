% Read data from files
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/RD09_finger2_doubleRes_subsampled;

clear all;
close all;

% Functions
norm_distance = @(x, y) sqrt(sum((x(:) - y(:)).*(x(:) - y(:))));

%========================================================================================================================
% NORMALIZE DATA
%========================================================================================================================
% Import data
filenameData = './input_data/forwardSignal_reference_3600sensors.dat';
timeSignal = importdata(filenameData, ' ', 0);
timeRD = timeSignal(1, :);
y0 = timeSignal(2:end, :);
% Plot
figure;
imagesc(y0);
box on;
colorbar();

% Normalize
NORM = 2.5e17;
y0_norm = y0/NORM;
timeSignal_norm = [timeRD; y0_norm];
% Plot
figure;
imagesc(y0_norm);
box on;
colorbar();

dlmwrite('./input_data/forwardSignal_reference_norm_3600sensors.dat', timeSignal_norm, 'delimiter', ' ');


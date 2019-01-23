% Read data from files
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/RD04_palm;

clear all;
close all;

% Functions
norm_distance = @(x, y) sqrt(sum((x(:) - y(:)).*(x(:) - y(:))));

%==================================================
% Dimensions
%==================================================
% Import dimensions
dim = importdata('input_data/dimensions.dat', ' ', 0);
Nx = dim(1, 1); dx = dim(2, 1);
Ny = dim(1, 2); dy = dim(2, 2);
Nz = dim(1, 3); dz = dim(2, 3);

%========================================================================================================================
% NORMALIZE DATA
%========================================================================================================================
% Import data
filenameData = './input_data/forwardSignal_reference_8736sensors.dat';
timeSignal = importdata(filenameData, ' ', 0);
timeRD = timeSignal(1, :);
y0 = timeSignal(2:end, :);
% Plot
figure;
imagesc(y0);
box on;
colorbar();

% Normalize
NORM = 4e16;
y0_norm = y0/NORM;
timeSignal_norm = [timeRD; y0_norm];
% Plot
figure;
imagesc(y0_norm);
box on;
colorbar();

dlmwrite('./input_data/forwardSignal_reference_norm_8736sensors.dat', timeSignal_norm, 'delimiter', ' ');

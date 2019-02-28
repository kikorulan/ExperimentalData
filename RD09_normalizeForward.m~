% Read data from files
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/RD08_finger2_doubleRes;

clear all;
close all;

% Functions
norm_distance = @(x, y) sqrt(sum((x(:) - y(:)).*(x(:) - y(:))));

%========================================================================================================================
% NORMALIZE DATA
%========================================================================================================================
% Import data
filenameData = './input_data/forwardSignal_reference_14400sensors.dat';
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

dlmwrite('./input_data/forwardSignal_reference_norm_14400sensors.dat', timeSignal_norm, 'delimiter', ' ');


%========================================================================================================================
% ANALIZE DATA
%========================================================================================================================

% Compute correlation coefficient
nrows_x = 120;
ncols_y = 120;
corr_coef = zeros(nrows_x-2, ncols_y-2);
for jj = 1:ncols_y-2
    for ii = 1:nrows_x-2
        index = ii + jj*nrows_x + 1;
        index_xpos = index + 1;
        index_xneg = index - 1;
        index_ypos = index + nrows_x;
        index_yneg = index - nrows_x;
        matcoef1 = corrcoef(y0_norm(index, :), y0_norm(index_xpos, :));
        matcoef2 = corrcoef(y0_norm(index, :), y0_norm(index_xneg, :));
        matcoef3 = corrcoef(y0_norm(index, :), y0_norm(index_ypos, :));
        matcoef4 = corrcoef(y0_norm(index, :), y0_norm(index_yneg, :));
        corr_coef(ii,jj) = (matcoef1(1, 2) + matcoef2(1, 2) + matcoef3(1, 2) + matcoef4(1, 2))/4;
    end
end
figure;
imagesc(corr_coef);

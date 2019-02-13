% Read data from files
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/RD06_finger2;

clear all;
close all;

% Functions
norm_distance = @(x, y) sqrt(sum((x(:) - y(:)).*(x(:) - y(:))));

%========================================================================================================================
% NORMALIZE DATA
%========================================================================================================================
% Import data
filenameData = './input_data/forwardSignal_reference_19152sensors.dat';
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

%dlmwrite('./input_data/forwardSignal_reference_norm_19152sensors.dat', timeSignal_norm, 'delimiter', ' ');



%========================================================================================================================
% ANALIZE DATA
%========================================================================================================================
nrows_x = 144;
ncols_y = 133;

% Mean and var
sensor_data = 2.5e17*y0_norm;
meanData = sum(sensor_data,2)/200;
varData = sum((sensor_data - repmat(meanData, [1 200])).^2,2)/200;
figure, imagesc(reshape(meanData, [144,133]))
figure, imagesc(reshape(varData, [144,133]))

% Compute correlation coefficient
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
        corr_coef(ii+1,jj+1) = (matcoef1(1, 2) + matcoef2(1, 2) + matcoef3(1, 2) + matcoef4(1, 2))/4;
    end
end

% Correlation coefficient
figure;
imagesc(corr_coef);

% Threshold 0.3
corr_coef_03 = corr_coef;
corr_coef_03(corr_coef_03 > 0.3) = 1;
corr_coef_03(corr_coef_03 <= 0.3) = 0;
figure;
imagesc(corr_coef_03);
colorbar();
% Percentage of valid sensors
perc = sum(corr_coef_03(:))/size(corr_coef_03(:), 1)


%=================================================================
% REDUCE DIMENSIONS
%=================================================================
% array_x = 1:144 -> array_x = 11:130
% array_y = 1:133 -> array_y =  3:123
x_min = 11;
x_max = 130;
y_min = 3;
y_max = 122;
corr_coef_03_reduce = corr_coef_03(x_min+1:x_max-1, y_min+1:y_max-1);
figure;
imagesc(corr_coef_03_reduce);
colorbar();
% Percentage of valid sensors
perc = sum(corr_coef_03_reduce(:))/size(corr_coef_03_reduce(:), 1)


% EVALUATE HOW MUCH DATA HAS BEEN DELETED
corr_coef_03_del = corr_coef_03;
corr_coef_03_del(x_min+1:x_max-1, y_min+1:y_max-1) = 0;
perc_red = sum(corr_coef_03_del(:))/size(corr_coef_03_del(:), 1)





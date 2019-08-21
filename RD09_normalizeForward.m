% Read data from files
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/RD09_finger2_doubleRes_subsampled;

clear all;
close all;

% Functions
norm_distance = @(x, y) sqrt(sum((x(:) - y(:)).*(x(:) - y(:))));

load ./settings/timeseriesdata_c1582_offset10_fulldata;
load ./settings/Full_scan1_temp@850nm_t0[0]_dx[106µm]_dy[106µm]_dt[17ns]_34s19m15h_10-05-18_avg1_2D_raw(QRS);



%========================================================================================================================
% NORMALIZE DATA
%========================================================================================================================
% Import data
filenameData = './input_data/forwardSignal_reference_14400sensors.dat';
%filenameData = './input_data/forwardSignal_reference_3600sensors.dat';
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
%dlmwrite('./input_data/forwardSignal_reference_norm_3600sensors.dat', timeSignal_norm, 'delimiter', ' ');

%========================================================================================================================
% CORRELATION COEFFICIENT
%========================================================================================================================
% Compute correlation coefficient
nrows_x = 144;
ncols_y = 133;
corr_coef = zeros(nrows_x-2, ncols_y-2);
for jj = 1:ncols_y-2
    for ii = 1:nrows_x-2
        ix = ii+1;
        iy = jj+1;
        ix_pos = ix+1;
        ix_neg = ix-1;
        iy_pos = iy+1;
        iy_neg = iy-1;
        matcoef1 = corrcoef(sensor_data(ix, iy, :), sensor_data(ix_pos, iy, :));
        matcoef2 = corrcoef(sensor_data(ix, iy, :), sensor_data(ix_neg, iy, :));
        matcoef3 = corrcoef(sensor_data(ix, iy, :), sensor_data(ix, iy_pos, :));
        matcoef4 = corrcoef(sensor_data(ix, iy, :), sensor_data(ix, iy_neg, :));
        corr_coef(ii,jj) = (matcoef1(1, 2) + matcoef2(1, 2) + matcoef3(1, 2) + matcoef4(1, 2))/4;
    end
end
figure;
imagesc(corr_coef);


%========================================================================================================================
% RMS
%========================================================================================================================
nrows_x = 144;
ncols_y = 133;

% RMS coeff with length 50
sizeVector = 50;
rms_coeff_50 = zeros(nrows_x, ncols_y);
for jj = 1:ncols_y
    for ii = 1:nrows_x
        sensor_index = sensor_data(ii, jj, end-sizeVector:end);
        rms = sum(sensor_index.*sensor_index);
        rms_coeff_50(ii, jj) = rms;
    end
end
figure;
imagesc(log(sqrt(rms_coeff_50)));
save ./input_data/rms_coeff_50 rms_coeff_50;

% RMS coeff with all steps
sizeVector = 489;
rms_coeff = zeros(nrows_x, ncols_y);
for jj = 1:ncols_y
    for ii = 1:nrows_x
        sensor_index = sensor_data(ii, jj, end-sizeVector:end);
        rms = sum(sensor_index.*sensor_index);
        rms_coeff(ii, jj) = rms;
    end
end
figure;
imagesc(log(sqrt(rms_coeff./rms_coeff_50)));

%=========================================================================
% PLOT
%=========================================================================
% ADJOINT - GENERATED PARAMETERS
PML_size = 10;
p0_adj_PML = h5read('./output_data/RD09_adjoint_output_14400sensors.h5', '/p_final');
p0_recon_adj = p0_adj_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
plot_projection(p0_recon_adj, dx);


%=========================================================================
% PLOT SLICE
%=========================================================================
for ii = 1:30
    p0_slice = reshape(p0_recon_adj(ii, :, :), [Ny, Nz]);
    p0_slice = p0_slice.*p0_slice;
    figure;
    imagesc(p0_slice);
end

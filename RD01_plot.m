% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/RD01_finger15;

clear all;
close all;

% Import dimensions
dim = importdata('input_data/dimensions.dat', ' ', 0);
Nx = dim(1, 1); dx = dim(2, 1); Ny = dim(1, 2); dy = dim(2, 2); Nz = dim(1, 3); dz = dim(2, 3);


%=========================================================================
% Forward data
%=========================================================================
load input_data/simulation_postproc;
sensor_data = importdata('input_data/forwardSignal_reference_14400sensors.dat', ' ', 0);

sensor_data_iter = h5read('output_data/RD01_forward_output_14400sensors.h5', '/p');
figure;
imagesc(sensor_data);
figure;
imagesc(8*sensor_data_iter);


%=========================================================================
% PLOT
%=========================================================================
set(0,'DefaultFigurePaperPositionMode','auto');
%============================
% KWAVE TR
%============================
PML_size = 10;
load output_data/p0_original;
p0_tr = p0_original(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
p0_recon_tr = max(0, p0_tr);
plot_projection(p0_recon_tr, dx);
saveas(gcf, 'output_data/RD01_kWave_tr', 'png');
saveas(gcf, 'output_data/RD01_kWave_tr.fig');

%============================
% KWAVE ADJOINT
%============================
PML_size = 10;
p0_adj_PML = h5read('output_data/RD01_adjoint_output_57600sensors.h5', '/p_final');
p0_adj = p0_adj_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
p0_recon_adj = max(0, p0_adj);
plot_projection(p0_recon_adj, dx);
saveas(gcf, 'output_data/RD01_kWave_adjoint', 'png');
saveas(gcf, 'output_data/RD01_kWave_adjoint.fig');

%============================
% RT ADJOINT
%============================
pixelPressureMatrix = importdata('output_data/PixelPressure.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
saveas(gcf, 'output_data/RD01_RT_adjoint', 'png');
saveas(gcf, 'output_data/RD01_RT_adjoint.fig');

%=========================================================================
% PLOT - SUBSAMPLED
%=========================================================================
%============================
% KWAVE ADJOINT - SUBSAMPLED
%============================
PML_size = 10;
p0_sub_PML = h5read('output_data/RD01_subsampled_output.h5', '/p_final');
p0_sub = p0_sub_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
p0_recon_sub = max(0, p0_sub);
plot_projection(p0_recon_sub);

%============================
% RT ADJOINT - SUBSAMPLED
%============================
pixelPressureMatrix = importdata('output_data/PixelPressure_subsampled.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, 240));
plot_projection(pixelPressure);


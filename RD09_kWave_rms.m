% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/RD09_finger2_doubleRes_subsampled;

clear all;
close all;

load ./settings/timeseriesdata_c1582_offset10;
load ./settings/Full_scan1_temp@850nm_t0[0]_dx[106µm]_dy[106µm]_dt[17ns]_34s19m15h_10-05-18_avg1_2D_raw(QRS);
%=========================================================================
% DOMAIN DEFINITION
%=========================================================================
% create the computational grid
Nx = 80;           % number of grid points in the x (row) direction
Ny = 288;           % number of grid points in the y (column) direction
Nz = 266;           % number of grid points in the y (column) direction
dx = 5.3e-5;        % grid point spacing in the x direction [m]
dy = 5.3e-5;        % grid point spacing in the y direction [m]
dz = 5.3e-5;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy, Nz, dz);

%==============================
% define the properties of the propagation medium
%==============================
% Load sound speed
c0 = 1582.01;
medium.sound_speed = c0;
medium.density = 1;

% compute time
dt = 1.6667e-8;
Nt = 200;
tMax = dt*(Nt-1);
kgrid.t_array = 0:dt:tMax;

%==============================
% EXPORT DATA TO TXT
%==============================
% REDUCE DIMENSIONS
% array_y = 1:144 -> array_y = 11:130
% array_z = 1:133 -> array_z =  3:123
y_min = 11;
y_max = 130;
z_min = 3;
z_max = 122;
% Reshape sensor data
load input_data/rms_coeff_50;
sensor_data = sensor_data./repmat(sqrt(rms_coeff_50), [1, 1, Nt]);
% Reshape sensor data
sensor_data = reshape(sensor_data, [Ny*Nz/2/2, Nt]);
forward_signal = [kgrid.t_array; sensor_data];
dlmwrite('input_data/forwardSignal_reference_14400sensors_rmsCorrection.dat', forward_signal, 'delimiter', ' ');

%=========================================================================
% SIMULATION
%=========================================================================
setenv LD_LIBRARY_PATH '/cs/research/medim/projects2/projects/frullan/lib/root/lib64:/cs/research/medim/projects2/projects/frullan/lib/glibc-2.27/mathvec';

%==============================
% ADJOINT - GENERATED PARAMETERS
%==============================
% set the input arguments: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};

% Consider all sensors
source_adj.p_mask = zeros(Nx, Ny, Nz);
source_adj.p_mask(1, 1:2:end, 1:2:end) = 1;
source_adj.p = fliplr(sensor_data);
%source_adj.p_mode = 'additive';
source_adj.p_mode = 'dirichlet';
% Sensor
sensor_adj.mask = ones(Nx, Ny, Nz);
sensor_adj.record = {'p_final'};
% Save and run
kspaceFirstOrder3D(kgrid, medium, source_adj, sensor_adj, input_args{:}, 'SaveToDisk', 'input_data/RD09_adjoint_input_14400sensors_rmsCorrection.h5');
system('../kspaceFirstOrder3D-CUDA -i input_data/RD09_adjoint_input_14400sensors_rmsCorrection.h5 -o output_data/RD09_adjoint_output_14400sensors_rmsCorrection.h5 --p_final');

%==============================
% FORWARD
%==============================
PML_size = 10;
p0_adj_PML = h5read('./output_data/RD09_adjoint_output_14400sensors_rmsCorrection.h5', '/p_final');
p0_adj = p0_adj_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);

% Source
source.p0 = smooth(kgrid, p0_adj, true);
source.p0 = max(0, source.p0);

% Define the sensors
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(1, 1:2:end, 1:2:end) = 1;
% Number of sensors
numberSensors = sum(sensor.mask(:))
save input_data/kgrid_data_3600sensors.mat kgrid medium source sensor input_args;
filename = './input_data/RD09_forward_input_14400sensors_rmsCorrection.h5';
kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', filename);
% Run forward
system('../kspaceFirstOrder3D-CUDA -i input_data/RD09_forward_input_14400sensors_rmsCorrection.h5 -o output_data/RD09_forward_output_14400sensors_rmsCorrection.h5');

%=========================================================================
% PLOT
%=========================================================================
% ADJOINT
PML_size = 10;
p0_adj_PML = h5read('./output_data/RD09_adjoint_output_14400sensors_rmsCorrection.h5', '/p_final');
p0_adj = p0_adj_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
p0_recon_adj = max(0, p0_adj);
plot_projection(p0_recon_adj, dx);

% ADJOINT - RMS CORRECTION
PML_size = 10;
p0_adj_PML = h5read('./output_data/RD09_adjoint_output_14400sensors.h5', '/p_final');
p0_adj = p0_adj_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
p0_recon_adj = max(0, p0_adj);
plot_projection(p0_recon_adj, dx);



% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/RD04_palm;

clear all;
close all;

load ./input_data/sensor_data_8736;

%=========================================================================
% DOMAIN DEFINITION
%=========================================================================
% create the computational grid
Nx = 50;           % number of grid points in the x (row) direction
Ny = 192;           % number of grid points in the y (column) direction
Nz = 182;           % number of grid points in the y (column) direction
dx = 5.5e-5;        % grid point spacing in the x direction [m]
dy = 5.5e-5;        % grid point spacing in the y direction [m]
dz = 5.5e-5;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy, Nz, dz);

%==============================
% define the properties of the propagation medium
%==============================
% Load sound speed
c0 = 1540.00001;
medium.sound_speed = c0;
medium.density = 1;

% compute time
dt = 1.7e-8;
Nt = 207;
tMax = dt*(Nt-1);
kgrid.t_array = 0:dt:tMax;

%==============================
% EXPORT DATA TO TXT
%==============================
pixelPressure = zeros(Nx*Nz, Ny);
dlmwrite('input_data/pixelPressure_0.dat', pixelPressure, 'delimiter', ' ');
sound_speed = c0*ones(Nx*Nz, Ny);
dlmwrite('input_data/sound_speed.dat', sound_speed, 'delimiter', ' ');
forward_signal = [kgrid.t_array; sensor_data];
dlmwrite('input_data/forwardSignal_reference_8736sensors.dat', forward_signal, 'delimiter', ' ');

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
source_adj.p_mode = 'additive';
 % Sensor
sensor_adj.mask = ones(Nx, Ny, Nz);
sensor_adj.record = {'p_final'};
% Save and run
kspaceFirstOrder3D(kgrid, medium, source_adj, sensor_adj, input_args{:}, 'SaveToDisk', 'input_data/RD04_adjoint_input_8736sensors.h5');
system('../kspaceFirstOrder3D-CUDA -i input_data/RD04_adjoint_input_8736sensors.h5 -o output_data/RD04_adjoint_output_8736sensors.h5 --p_final');

%==============================
% FORWARD
%==============================
PML_size = 10;
p0_adj_PML = h5read('output_data/RD04_adjoint_output_8736sensors.h5', '/p_final');
p0_adj = p0_adj_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
% Source
source.p0 = smooth(kgrid, p0_adj, true);
source.p0 = max(0, source.p0);
% Define the sensors
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(1, 1:2:end, 1:2:end) = 1;
% Number of sensors
numberSensors = sum(sensor.mask(:))
save input_data/kgrid_data_8736sensors.mat kgrid medium source sensor input_args;
filename = 'input_data/RD04_forward_input_8736sensors.h5';
kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', filename);
% Run forward
system('../kspaceFirstOrder3D-OMP -i input_data/RD04_forward_input_8736sensors.h5 -o output_data/RD04_forward_output_8736sensors.h5');


%=========================================================================
% SENSOR DATA
%=========================================================================
sensor_data_forward_8736 = h5read('output_data/RD04_forward_output_8736sensors.h5', '/p');

%=========================================================================
% PLOT
%=========================================================================
% ADJOINT - GENERATED PARAMETERS
PML_size = 10;
p0_adj_PML = h5read('output_data/RD04_adjoint_output_8736sensors.h5', '/p_final');
p0_adj = p0_adj_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
p0_recon_adj = max(0, p0_adj);
plot_projection(p0_recon_adj, dx);



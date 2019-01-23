% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/RD01_finger15;

clear all;
close all;

cd input_data;
load simulation_postproc;
cd ..;

%=========================================================================
% DOMAIN DEFINITION
%=========================================================================
% create the computational grid
Nx = 80;           % number of grid points in the x (row) direction
Ny = 240;           % number of grid points in the y (column) direction
Nz = 240;           % number of grid points in the y (column) direction
dx = 5.3e-5;        % grid point spacing in the x direction [m]
dy = 5.3e-5;        % grid point spacing in the y direction [m]
dz = 5.3e-5;        % grid point spacing in the y direction [m]
kgrid_adj = makeGrid(Nx, dx, Ny, dy, Nz, dz);

%==============================
% define the properties of the propagation medium
%==============================
% Load sound speed
c0 = 1580.00001;
medium_adj.sound_speed = c0;
medium_adj.density = 1;

% compute time
dt = 1.6667e-8;
Nt = 486;
tMax = dt*(Nt-1);
kgrid_adj.t_array = 0:dt:tMax;

%==============================
% EXPORT DATA TO TXT
%==============================
sound_speed = c0*ones(Nx*Nz, Ny);
dlmwrite('input_data/sound_speed.dat', sound_speed, 'delimiter', ' ');
forward_signal = [kgrid.t_array; sensor_data_postproc.time_reversal_boundary_data];
dlmwrite('input_data/forwardSignal_reference_57600sensors.dat', forward_signal, 'delimiter', ' ');


%=========================================================================
% SIMULATION
%=========================================================================
setenv LD_LIBRARY_PATH '/cs/research/medim/projects2/projects/frullan/lib/root/lib64';

%==============================
% ADJOINT - GENERATED PARAMETERS
%==============================
% set the input arguments: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};

% Consider all sensors
source_adj.p_mask = zeros(Nx, Ny, Nz);
source_adj.p_mask(1, :, :) = 1;
source_adj.p = fliplr(sensor_data_postproc.time_reversal_boundary_data);
source_adj.p_mode = 'additive';
 % Sensor
sensor_adj.mask = ones(Nx, Ny, Nz);
sensor_adj.record = {'p_final'};
% Save and run
kspaceFirstOrder3D(kgrid_adj, medium_adj, source_adj, sensor_adj, input_args{:}, 'SaveToDisk', 'input_data/RD01_adjoint_input_57600sensors.h5');
system('../kspaceFirstOrder3D-OMP -i input_data/RD01_adjoint_input_57600sensors.h5 -o output_data/RD01_adjoint_output_57600sensors.h5 --p_final');

%=========================================================================
% PLOT
%=========================================================================

% ADJOINT - GENERATED PARAMETERS
PML_size = 10;
p0_adj_PML = h5read('output_data/RD01_adjoint_output_57600sensors.h5', '/p_final');
p0_adj = p0_adj_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
p0_recon_adj = max(0, p0_adj);
plot_projection(p0_recon_adj);



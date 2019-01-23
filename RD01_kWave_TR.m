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
kgrid_tr = makeGrid(Nx, dx, Ny, dy, Nz, dz);

%==============================
% define the properties of the propagation medium
%==============================
% Load sound speed
c0 = 1580.00001;
medium_tr.sound_speed = c0;
medium_tr.density = 1;

% compute time
dt = 1.6667e-8;
Nt = 486;
tMax = dt*(Nt-1);
kgrid_tr.t_array = 0:dt:tMax;


%=========================================================================
% SIMULATION
%=========================================================================
setenv LD_LIBRARY_PATH '/cs/research/medim/projects2/projects/frullan/lib/root/lib64';
%==============================
% TIME REVERSAL - ORIGINAL PARAMETERS
%==============================
p0_original = kspaceFirstOrder3DC(kgrid, medium, source, sensor_data_postproc, input_args_here{:});
save output_data/p0_original p0_original;

%==============================
% TIME REVERSAL - Compute forward pressure
%==============================
% set the input arguments: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};

% Build initial pressure
clear source;
source.p0 = p0_original;
 % Sensor
clear sensor;
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(1, :, :) = 1;

% Save and run
kspaceFirstOrder3D(kgrid_tr, medium_tr, source, sensor, input_args{:}, 'SaveToDisk', 'input_data/RD01_tr_input.h5');
system('../kspaceFirstOrder3D-OMP -i input_data/RD01_tr_input.h5 -o output_data/RD01_tr_output.h5 --p_final');

%=========================================================================
% PLOT
%=========================================================================

% TR p0 original
PML_size = 10;
load ./output_data/p0_original;
p0 = p0_original(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
p0_recon = max(0, p0);
plot_projection(p0_recon, dx);





% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/RD01_finger15;

clear all;
close all;

load input_data/simulation_postproc;

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
kgrid_sub = makeGrid(Nx, dx, Ny, dy, Nz, dz);

%==============================
% define the properties of the propagation medium
%==============================
% Load sound speed
c0 = 1580.00001;
medium.sound_speed = c0;
medium.density = 1;

% compute time
dt = 1.6667e-8;
Nt = 486;
tMax = dt*(Nt-1);
kgrid.t_array = 0:dt:tMax;

% set the input arguments: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};
%==============================
% Forward
%==============================
nSensorsArray = 120;
% Define the sensors
sensor.mask = zeros(Nx, Ny, Nz);
xArray = round(1:(Nx-1)/(nSensorsArray-1):Nx);
yArray = round(1:(Ny-1)/(nSensorsArray-1):Ny);
zArray = round(1:(Nz-1)/(nSensorsArray-1):Nz);
% YZ faces
[Y, Z] = meshgrid(yArray, zArray);
sensorYZ1 = 1    + Nx*(Y(:)-1) + Nx*Ny*(Z(:)-1);
sensor.mask(sensorYZ1) = 1;
% Number of sensors
numberSensors = sum(sensor.mask(:))
% Linear array
factor = floor(Ny/nSensorsArray);
linIndex = 1:Ny*Nz;
linIndex_mat = reshape(linIndex, [Ny, Nz]);
linIndex_sub = linIndex_mat(1:factor:end, 1:factor:end);
sensor_data = sensor_data_postproc.time_reversal_boundary_data(linIndex_sub(:), :);
save input_data/sensor_data_14400sensors sensor_data;
save input_data/kgrid_data_14400sensors.mat kgrid medium sensor input_args;
%==============================
% EXPORT DATA TO TXT
%==============================
forward_signal = [kgrid.t_array; sensor_data];
dlmwrite('input_data/forwardSignal_reference_14400sensors.dat', forward_signal, 'delimiter', ' ');

%=========================================================================
% SIMULATION - ADJOINT
%=========================================================================
setenv LD_LIBRARY_PATH '/cs/research/medim/projects2/projects/frullan/lib/root/lib64:/cs/research/medim/projects2/projects/frullan/lib/glibc-2.27/mathvec';

D = gpuDevice(2)
% Consider all sensors
source.p_mask = zeros(Nx, Ny, Nz);
source.p_mask(1, 1:factor:end, 1:factor:end) = 1;
source.p = fliplr(sensor_data);
source.p_mode = 'additive';
 % Sensor
sensor.mask = ones(Nx, Ny, Nz);
sensor.record = {'p_final'};
% Save and run
kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', 'input_data/RD01_adjoint_input_14400sensors.h5');
system('../kspaceFirstOrder3D-CUDA -i input_data/RD01_adjoint_input_14400sensors.h5 -o output_data/RD01_adjoint_output_14400sensors.h5 --p_final -g 1');


%=========================================================================
% SIMULATION - FORWARD
%=========================================================================
% ADJOINT - 14400 sensors
PML_size = 10;
p0_PML = h5read('output_data/RD01_adjoint_output_14400sensors.h5', '/p_final');
p0 = max(0, p0_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size));
% Source
load ./input_data/kgrid_data_14400sensors;
clear source;
source.p0 = smooth(kgrid, p0, true);
source.p0 = max(0, source.p0);
% Save to disk
filename = 'input_data/RD01_forward_input_14400sensors.h5';
kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', filename);
% Call C++ code
system('../kspaceFirstOrder3D-CUDA -i input_data/RD01_forward_input_14400sensors.h5 -o output_data/RD01_forward_output_14400sensors.h5 -g 1');

%=========================================================================
% PLOT
%=========================================================================
% ADJOINT - 14400 sensors
PML_size = 10;
p0_sub_PML = h5read('output_data/RD01_adjoint_output_14400sensors.h5', '/p_final');
p0_sub = p0_sub_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
p0_recon_sub = max(0, p0_sub);
plot_projection(p0_recon_sub, dx);
% ADJOINT - 57600 sensors
PML_size = 10;
p0_sub_PML = h5read('output_data/RD01_adjoint_output_57600sensors.h5', '/p_final');
p0_sub = p0_sub_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
p0_recon_sub = max(0, p0_sub);
plot_projection(p0_recon_sub, dx);

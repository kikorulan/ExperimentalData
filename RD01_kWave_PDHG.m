% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/RD01_finger15;

clear all;
close all;

% LIBRARIES
setenv LD_LIBRARY_PATH '/cs/research/medim/projects2/projects/frullan/lib/root/lib64:/cs/research/medim/projects2/projects/frullan/lib/glibc-2.27/mathvec';

%====================
% LOAD DATA
%====================
load ./input_data/kgrid_data_14400sensors.mat;
load ./input_data/sensor_data_14400sensors;
y_0 = sensor_data;
% smooth the initial pressure distribution and restore the magnitude
PML_size = 10;
u_PML = h5read('output_data/RD01_adjoint_output_57600sensors.h5', '/p_final');
xBar_k = max(0, u_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size));
x_k   = xBar_k;
x_k_1 = xBar_k;
y_k   = y_0;
y_k_1 = y_0;

%============================================================
% SIMULATION
%============================================================

% PARAMETERS
tau = 1;
sigma = 0.5;
theta = 1;
lambda = 1e-4;
nIter = 10;

para.maxIter = 200;

for n = 1:nIter
    disp(n)
    %====================
    % FORWARD
    %====================
    clear source;
    source.p0 = smooth(kgrid, xBar_k, true);
    source.p0 = max(0, source.p0);
    % Save to disk
    filename = 'input_data/RD01_forward_input_A.h5';
    kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', filename);
    % Call C++ code
    system('../kspaceFirstOrder3D-CUDA -i input_data/RD01_forward_input_A.h5 -o output_data/RD01_forward_output_A.h5 -g 1');
    
    %==============================
    % ADJOINT
    %==============================
    % Read results
    ax_k = h5read('output_data/RD01_forward_output_A.h5', '/p');
    y_k = (y_k_1 + sigma*ax_k - sigma*y_0)/(1+sigma);
    % Number of sensors
    sensor_index = find(sensor.mask == 1);
    nSensors = length(sensor_index);
    
    % Consider all sensors
    source_adjoint.p_mask = sensor.mask;
    source_adjoint.p = fliplr(y_k);
    % Sensor
    sensor_adjoint.mask = ones(kgrid.Nx, kgrid.Ny, kgrid.Nz);
    sensor_adjoint.record = {'p_final'};
    % Save and run
    kspaceFirstOrder3D(kgrid, medium, source_adjoint, sensor_adjoint, input_args{:}, 'SaveToDisk', 'input_data/RD01_adjoint_input_A.h5');
    system('../kspaceFirstOrder3D-CUDA -i input_data/RD01_adjoint_input_A.h5 -o output_data/RD01_adjoint_output_A.h5 --p_final -g 1');
    
    %==============================
    % UPDATE 
    %==============================
    p_PML = h5read('output_data/RD01_adjoint_output_A.h5', '/p_final');
    p = p_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
    % x_k
    x_k = conTVdenoising(x_k - 2*tau*p, lambda, para);
    %x_k = softThresh(x_k - tau*p, lambda);
    x_k = max(0, x_k);
    xBar_k = x_k + theta*(x_k - x_k_1);
    % Update
    y_k_1 = y_k;
    x_k_1 = x_k;
    
    %==============================
    % PLOT
    %==============================
    plot_projection(x_k, kgrid.dx);
    colorbar();
    %pause(0.1);
end


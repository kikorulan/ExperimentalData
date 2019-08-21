% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/RD09_finger2_doubleRes_subsampled;


clear all;
close all;


% LIBRARIES
setenv LD_LIBRARY_PATH '/cs/research/medim/projects2/projects/frullan/lib/root/lib64';

%====================
% LOAD DATA
%====================
load ./input_data/kgrid_data_3600sensors.mat;
load ./input_data/sensor_data_3600.mat;

u_k = zeros(kgrid.Nx, kgrid.Ny, kgrid.Nz);
%============================================================
% SIMULATION
%============================================================
PML_size = 10;

% PARAMETERS
tau = 10;
lambda = 1e-4;
nIter = 10;

para.maxIter = 200;

for n = 1:nIter
    %====================
    % FORWARD
    %====================
    clear source;
    source.p0 = smooth(kgrid, u_k, true);
    source.p0 = max(0, source.p0);
    % Save to disk
    filename = './input_data/RD09_forward_input_A.h5';
    kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', filename);
    % Call C++ code
    system('../kspaceFirstOrder3D-OMP -i input_data/RD09_forward_input_A.h5 -o output_data/RD09_forward_output_A.h5');
    
    %==============================
    % ADJOINT
    %==============================
    % Read results
    sensor_data_adjoint = h5read('./output_data/RD09_forward_output_A.h5', '/p');
    % Consider all sensors
    source_adj.p_mask = sensor.mask;
    source_adj.p = fliplr(sensor_data_adjoint - sensor_data);
    % Sensor
    sensor_adj.record = {'p_final'};
    % Save and run
    kspaceFirstOrder3D(kgrid, medium, source_adj, sensor_adj, input_args{:}, 'SaveToDisk', './input_data/RD09_adjoint_input_A.h5');
    system('../kspaceFirstOrder3D-OMP -i input_data/RD09_adjoint_input_A.h5 -o output_data/RD09_adjoint_output_A.h5 --p_final');
    
    %==============================
    % UPDATE 
    %==============================
    p_PML = h5read('./output_data/RD09_adjoint_output_A.h5', '/p_final');
    p = p_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
    u_k = conTVdenoising(u_k - 2*tau*p, lambda, para);
    u_k = max(0, u_k);

    %==============================
    % PLOT
    %==============================
    plot_projection(u_k, kgrid.dx);
    pause(0.1);
end


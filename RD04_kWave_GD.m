% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/RD04_palm;

clear all;
close all;

% LIBRARIES
setenv LD_LIBRARY_PATH '/cs/research/medim/projects2/projects/frullan/lib/root/lib64:/cs/research/medim/projects2/projects/frullan/lib/glibc-2.27/mathvec';
% Distance function
norm_distance = @(x, y) sqrt(sum((x(:) - y(:)).*(x(:) - y(:))));

%====================
% LOAD DATA
%====================
load ./input_data/kgrid_data_8736sensors;
load ./input_data/sensor_data_8736;
% smooth the initial pressure distribution and restore the magnitude
PML_size = 10;
u_PML = h5read('./output_data/RD04_adjoint_output_8736sensors.h5', '/p_final');
u_k_neg = 8*u_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
u_k = max(0, u_k_neg);
plot_projection(u_k, 1);

figure;
imagesc(sensor_data);
title('Sensor data (real)');
colorbar();
% Save pressure
%%  u_matrix = cube2matrix(u_k);
%%  dlmwrite('input_data/pressure_adjoint_kWave_8736sensors.dat', u_matrix, 'delimiter', ' ');


%============================================================
% SIMULATION
%============================================================
% PARAMETERS
tau = 3e-1;
lambda = 1e-5;
thresh = 1e-4;
nIter = 50;

para.maxIter = 200;

forward_error = [];
for n = 1:nIter
    disp(n)
    %====================
    % FORWARD
    %====================
    clear source;
    source.p0 = smooth(kgrid, u_k, true);
    source.p0 = max(0, source.p0);
    % Save to disk
    filename = './input_data/RD04_forward_input_A.h5';
    kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', filename);
    % Call C++ code
    system('../kspaceFirstOrder3D-CUDA -i input_data/RD04_forward_input_A.h5 -o output_data/RD04_forward_output_A.h5 -g 1');
    
    %==============================
    % ADJOINT
    %==============================
    % Read results
    sensor_data_adjoint = h5read('./output_data/RD04_forward_output_A.h5', '/p');

    % Number of sensors
    sensor_index = find(sensor.mask == 1);
    nSensors = length(sensor_index);
    
    % Consider all sensors
    source_adjoint.p_mask = sensor.mask;
    source_adjoint.p = fliplr(sensor_data_adjoint - sensor_data);
    forward_error = [forward_error norm_distance(sensor_data_adjoint, sensor_data)];
    % Sensor
    sensor_adjoint.mask = ones(kgrid.Nx, kgrid.Ny, kgrid.Nz);
    sensor_adjoint.record = {'p_final'};
    % Save and run
    kspaceFirstOrder3D(kgrid, medium, source_adjoint, sensor_adjoint, input_args{:}, 'SaveToDisk', 'input_data/RD04_adjoint_input_A.h5');
    system('../kspaceFirstOrder3D-CUDA -i input_data/RD04_adjoint_input_A.h5 -o output_data/RD04_adjoint_output_A.h5 --p_final -g 1');
    
    %==============================
    % UPDATE 
    %==============================
    p_PML = h5read('./output_data/RD04_adjoint_output_A.h5', '/p_final');
    p = p_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
    u_k = conTVdenoising(u_k - 2*tau*p, lambda, para);
    %u_k = softThresh(u_k - 2*tau*p, thresh);
    u_k = max(0, u_k);

    %==============================
    % PLOT
    %==============================
    if (mod(n, 5) == 0)
        % FORWARD
        figure;
        imagesc(sensor_data_adjoint);
        title(['Sensor data iter ', int2str(n)]);
        colorbar();
        % ADJOINT
        plot_projection(u_k, kgrid.dx);
        title(['Projection iter ', int2str(n)]);
        colorbar();
        % Error
        figure;
        plot(forward_error);
    end
    %pause(0.1);
end



%============================================================
% SAVE
%============================================================
u_matrix = cube2matrix(u_k);
dlmwrite('output_data/pressure_kWave_GD_200iter.dat', u_matrix, 'delimiter', ' ');

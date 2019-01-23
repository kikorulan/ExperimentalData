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
u_k_neg = u_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
y_k = max(0, u_k_neg);
u_k_1 = y_k;
plot_projection(y_k, 1);

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
t_k = 1;
t_k_1 = 1;
tau = 1e-1;
lambda = 5e-5;
nIter = 200;

para.maxIter = 200;

forward_error = [];
for n = 1:nIter
    disp(n)
    %====================
    % FORWARD
    %====================
    clear source;
    source.p0 = smooth(kgrid, y_k, true);
    source.p0 = max(0, source.p0);
    % Save to disk
    filename = './input_data/RD04_forward_input_FISTA_A.h5';
    kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', filename);
    % Call C++ code
    system('../kspaceFirstOrder3D-CUDA -i input_data/RD04_forward_input_FISTA_A.h5 -o output_data/RD04_forward_output_FISTA_A.h5 -g 0');
    
    %==============================
    % ADJOINT
    %==============================
    % Read results
    sensor_data_adjoint = h5read('./output_data/RD04_forward_output_FISTA_A.h5', '/p');

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
    kspaceFirstOrder3D(kgrid, medium, source_adjoint, sensor_adjoint, input_args{:}, 'SaveToDisk', 'input_data/RD04_adjoint_input_FISTA_A.h5');
    system('../kspaceFirstOrder3D-CUDA -i input_data/RD04_adjoint_input_FISTA_A.h5 -o output_data/RD04_adjoint_output_FISTA_A.h5 --p_final -g 0');
    
    %==============================
    % UPDATE 
    %==============================
    p_PML = h5read('./output_data/RD04_adjoint_output_FISTA_A.h5', '/p_final');
    p = p_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
    u_k = conTVdenoising(y_k - 2*tau*p, lambda, para);
    u_k = max(0, u_k);
    % t_k
    t_k = (1 + sqrt(1 + 4*t_k_1*t_k_1))/2;
    % y_k
    y_k = u_k + (t_k_1 - 1)/t_k*(u_k - u_k_1);
    % Next iter
    u_k_1 = u_k;
    t_k_1 = t_k;

    %==============================
    % PLOT
    %==============================
    if (mod(n, 20) == 0)
        % FORWARD
        figure;
        imagesc(sensor_data_adjoint);
        title(['(FISTA) Sensor data iter ', int2str(n)]);
        colorbar();
        % ADJOINT
        plot_projection(u_k, kgrid.dx);
        title(['(FISTA) Projection iter ', int2str(n)]);
        colorbar();
    end
end



%============================================================
% SAVE
%============================================================
u_matrix = cube2matrix(u_k);
dlmwrite('output_data/pressure_kWave_FISTA_200iter.dat', u_matrix, 'delimiter', ' ');

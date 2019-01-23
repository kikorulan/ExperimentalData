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
y_0 = sensor_data;
% smooth the initial pressure distribution and restore the magnitude
PML_size = 10;
u_PML = h5read('./output_data/RD04_adjoint_output_8736sensors.h5', '/p_final');
xBar_k= max(0, u_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size));
x_k   = xBar_k;
x_k_1 = xBar_k;
y_k   = y_0;
y_k_1 = y_0;


%============================================================
% SIMULATION
%============================================================
% PARAMETERS
tau = 1e-1;
sigma = 1;
theta = 1;
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
    source.p0 = smooth(kgrid, xBar_k, true);
    source.p0 = max(0, source.p0);
    % Save to disk
    filename = './input_data/RD04_forward_input_PDHG_A.h5';
    kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', filename);
    % Call C++ code
    system('../kspaceFirstOrder3D-CUDA -i input_data/RD04_forward_input_PDHG_A.h5 -o output_data/RD04_forward_output_PDHG_A.h5 -g 0');
    
    %==============================
    % ADJOINT
    %==============================
    % Read results
    ax_k = h5read('./output_data/RD04_forward_output_PDHG_A.h5', '/p');
    y_k = (y_k_1 + sigma*ax_k - sigma*y_0)/(1+sigma);
    % Number of sensors
    sensor_index = find(sensor.mask == 1);
    nSensors = length(sensor_index);
    
    % Consider all sensors
    source_adjoint.p_mask = sensor.mask;
    source_adjoint.p = fliplr(y_k);
    forward_error = [forward_error norm_distance(ax_k, sensor_data)];
    % Sensor
    sensor_adjoint.mask = ones(kgrid.Nx, kgrid.Ny, kgrid.Nz);
    sensor_adjoint.record = {'p_final'};
    % Save and run
    kspaceFirstOrder3D(kgrid, medium, source_adjoint, sensor_adjoint, input_args{:}, 'SaveToDisk', 'input_data/RD04_adjoint_input_PDHG_A.h5');
    system('../kspaceFirstOrder3D-CUDA -i input_data/RD04_adjoint_input_PDHG_A.h5 -o output_data/RD04_adjoint_output_PDHG_A.h5 --p_final -g 0');
    
    %==============================
    % UPDATE 
    %==============================
    p_PML = h5read('./output_data/RD04_adjoint_output_PDHG_A.h5', '/p_final');
    p = p_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
    % x_k
    x_k = conTVdenoising(x_k - tau*p, lambda, para);
    x_k = max(0, x_k);
    xBar_k = x_k + theta*(x_k - x_k_1);
    % Update
    y_k_1 = y_k;
    x_k_1 = x_k;

    %==============================
    % PLOT
    %==============================
    if (mod(n, 20) == 0)
        % FORWARD
        figure;
        imagesc(ax_k);
        title(['(PDHG) Sensor data iter ', int2str(n)]);
        colorbar();
        % ADJOINT
        plot_projection(x_k, kgrid.dx);
        title(['(PDHG) Projection iter ', int2str(n)]);
        colorbar();
    end
end



%============================================================
% SAVE
%============================================================
u_matrix = cube2matrix(u_k);
dlmwrite('output_data/pressure_kWave_FISTA_200iter.dat', u_matrix, 'delimiter', ' ');

% Read data from files
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/RD06_finger2;

clear all;
close all;


% Functions
norm_distance = @(x, y) sqrt(sum((x(:) - y(:)).*(x(:) - y(:))));

%==================================================
% Dimensions
%==================================================
% Import dimensions
dim = importdata('input_data/dimensions.dat', ' ', 0);
Nx = dim(1, 1); dx = dim(2, 1);
Ny = dim(1, 2); dy = dim(2, 2);
Nz = dim(1, 3); dz = dim(2, 3);

%============================================================================================================================================
% FORWARD PROBLEM
%============================================================================================================================================

%==================================================
% TIME SIGNAL - REAL DATA
%==================================================
% Import data
filenameData = './input_data/forwardSignal_reference_19152sensors.dat';
timeSignal = importdata(filenameData, ' ', 0);
timeRD = timeSignal(1, :);
y0_raw = timeSignal(2:end, :);
% Plot
figure;
imagesc(y0_raw);
box on;
colorbar();

% Import data
filenameData = './input_data/forwardSignal_reference_norm_19152sensors.dat';
timeSignal = importdata(filenameData, ' ', 0);
timeRD = timeSignal(1, :);
y0 = timeSignal(2:end, :);
% Plot
figure;
imagesc(y0);
box on;
colorbar();

%==================================================
% TIME SIGNAL - kWave
%==================================================
sensor_data = h5read('./output_data/RD04_forward_output_19152sensors.h5', '/p');
inputKWave = sensor_data(:, :);
% Plot
figure;
imagesc(inputKWave);
box on;
colorbar();

%============================================================================================================================================
% ADJOINT PROBLEM
%============================================================================================================================================

%==================================================
% ADJOINT k-Wave
%==================================================
PML_size = 10;
pixelPressureMatrix = h5read('./output_data/RD04_adjoint_output_19152sensors.h5', '/p_final');
pixelPressure = max(0, pixelPressureMatrix(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size));
plot_projection(pixelPressure, dx);

%==================================================
% ADJOINT RT
%==================================================
pixelPressureMatrix = importdata('output_data/pixelPressure_adjoint_RT.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

%========================================================================================================================
% FORWARD DATA
%========================================================================================================================
%==============================
% Stochastic Gradient Descent
%==============================
%SGD.tau = {'1e16', '2e16', '5e16', '1e17', '2e17', '5e17', '1e18', '2e18', '5e18'};
SGD.tau = {'5e17'};
SGD.lambda = {'1e-4', '2e-4', '4e-4', '8e-4'};
%SGD.lambda = {'1e-4', '2e-4', '4e-4', '8e-4'};

clear yi_18;
SGD.epoch  = int2str(1);
for j = 1:length(SGD.lambda)
    for i = 1:length(SGD.tau)
        tSignal = importdata(['./results/forward/S-FB/forwardSignal_S-GD_tau', SGD.tau{i}, '_lambda', SGD.lambda{j}, '_batch', SGD.batch, '_epoch', SGD.epoch, '.dat'], ' ', 0);
        yi_18{j} = tSignal(2:end, :);
        figure;
        imagesc(yi_18{j})
        t = title(['S-GD - t = ', SGD.tau{i}, ', l = ', SGD.lambda{j}, ', batch = ', SGD.batch, ', epoch = ', SGD.epoch, ' - homogeneous SS']);
        t.Visible = 'on'; 
        pause(0.1);
    end
end

%========================================================================================================================
% ITERATIVE RECONSTRUCTION
%========================================================================================================================

iter = 10;
%============================================================
% PARAMETERS
%============================================================
% GD *************************************
GD = [];
GD.tau    = '5e17';
GD.lambda = '2e-3';
GD.iter   = int2str(iter);
% S-GD ***********************************
SGD = [];
SGD.tau    = '5e17';
SGD.lambda = '1e-5';
SGD.batch  = '90';
SGD.epoch  = int2str(iter);
% FISTA **********************************
FISTA = [];
FISTA.tau    = '2e17';
FISTA.lambda = '2e-3';
FISTA.iter   = int2str(iter);
% PDHG ***********************************
PDHG = [];
PDHG.sigma  = '2e-1';
PDHG.tau    = '1e18';
PDHG.theta  = '1';
PDHG.lambda = '2e-3';
PDHG.iter   = int2str(iter);
% S-PDHG *********************************
SPDHG = [];
SPDHG.sigma  = '1';
SPDHG.tau    = '5e17';
SPDHG.theta  = '1';
SPDHG.lambda = '2e-5';
SPDHG.batch  = '90';
SPDHG.epoch  = int2str(iter);


%==============================
% Gradient Descent
%==============================
%GD.tau = {'1e16', '2e16', '5e16', '1e17', '2e17', '5e17', '1e18', '2e18'};
GD.tau = {'1e18'};
GD.lambda = {'2e-3', '5e-3', '1e-2'};
GD.iter = '5';
for j = 1:length(GD.lambda)
    for i = 1:length(GD.tau)
        pixelPressureMatrix = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau{i}, '_lambda', GD.lambda{j}, '_iter', GD.iter, '.dat'], ' ', 0);
        pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
        plot_projection(pixelPressure, dx);
        a = axes;
        t = title(['GD - tau = ', GD.tau{i}, ', lambda = ', GD.lambda{j}, ', iter = ', GD.iter]);
        a.Visible = 'off'; 
        t.Visible = 'on'; 
    end
    %saveas(gcf, ['./figures/RD06_GD_tau', GD.tau{i}, '_lambda', GD.lambda{j}, '_iter', GD.iter, '.fig']);
end

%==============================
% Stochastic Gradient Descent
%==============================
%SGD.tau = {'1e16', '2e16', '5e16', '1e17', '2e17', '5e17', '1e18', '2e18', '5e18'};
SGD.tau = {'2e18'};
SGD.lambda = {'1e-2', '2e-2'};
%SGD.lambda = {'5e-5', '1e-4', '2e-4', '5e-4', '1e-3'};
%SGD.lambda = {'1e-4', '2e-4', '4e-4', '8e-4'};
SGD.batch = '3800';

SGD.epoch  = int2str(5);
for j = 1:length(SGD.lambda)
    for i = 1:length(SGD.tau)
        pixelPressureMatrix = importdata(['./results/adjoint/S-FB/pixelPressure_S-GD_tau', SGD.tau{i}, '_lambda', SGD.lambda{j}, '_batch', SGD.batch, '_epoch', SGD.epoch, '.dat'], ' ', 0);
        pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
        plot_projection(pixelPressure, dx);
        a = axes;
        %colorbar();
        t = title(['S-GD - t = ', SGD.tau{i}, ', l = ', SGD.lambda{j}, ', batch = ', SGD.batch, ', epoch = ', SGD.epoch, ' - homogeneous SS']);
        a.Visible = 'off'; 
        t.Visible = 'on'; 
        pause(0.1);
    end
end
%saveas(gcf, ['./figures/Example74_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_epoch', SGD.epoch, '.fig']);

%==============================
% FISTA
%==============================
%%  pixelPressureMatrix = importdata(['./results/pixelPressure_FISTA_tau', FISTA.tau, '_lambda', FISTA.lambda, '_iter', FISTA.iter, '.dat'], ' ', 0);
%%  pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
%%  plot_projection(pixelPressure, dx);
%%  a = axes;
%%  t = title(['FISTA - t = ', FISTA.tau, ', l = ', FISTA.lambda, ', iter = ', FISTA.iter, ' - homogeneous SS']);
%%  a.Visible = 'off'; 
%%  t.Visible = 'on'; 
%%  %saveas(gcf, ['./figures/Example74_FISTA_tau', FISTA.tau, '_lambda', FISTA.lambda, '_iter', FISTA.iter, '.fig']);

%==============================
% PDHG
%==============================
PDHG.tau = {'5e17'};
PDHG.lambda = {'2e-3'};
PDHG.sigma = {'5e-1', '1'};
for ii = 1:length(PDHG.sigma)
    pixelPressureMatrix = importdata(['./results/adjoint/PDHG/pixelPressure_PDHG_sigma', PDHG.sigma{ii}, '_tau', PDHG.tau{1}, '_theta', PDHG.theta, '_lambda', PDHG.lambda{1}, '_iter', int2str(iter), '.dat'], ' ', 0);
    pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
    plot_projection(pixelPressure, dx);
    a = axes;
    t = title(['PDHG - s = ', PDHG.sigma{ii}, ', t = ', PDHG.tau{1}, ', l = ', PDHG.lambda{1}, ', iter = ', PDHG.iter, ' - homogeneous SS']);
    a.Visible = 'off'; 
    t.Visible = 'on'; 
end
%saveas(gcf, ['./figures/Example74_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau, '_lambda', PDHG.lambda, '_iter', PDHG.iter, '.fig']);

%==============================
% S-PDHG
%==============================
%%  pixelPressureMatrix = importdata(['./results/pixelPressure_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau, '_theta', SPDHG.theta, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_epoch', SPDHG.epoch, '.dat'], ' ', 0);
%%  pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
%%  plot_projection(pixelPressure, dx);
%%  a = axes;
%%  t = title(['S-PDHG - s = ', SPDHG.sigma, ', t = ', SPDHG.tau, ', l = ', SPDHG.lambda, ', batch = ', SPDHG.batch, ', epoch = ', SPDHG.epoch, ' - homogeneous SS']);
%%  a.Visible = 'off'; 
%%  t.Visible = 'on'; 
%%  %saveas(gcf, ['./figures/Example74_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_epoch', SPDHG.epoch, '.fig']);


%========================================================================================================================
% DUAL DISTANCE
%========================================================================================================================

%====================
% Gradient descent - tau variation
%====================
disp('GD');
GD.tau    = {'1e16', '2e16', '5e16', '1e17', '2e17', '5e17'};
GD.lambda = {'2e-3'};
clear GD_error_dd;
for l = 1:length(GD.tau)
    disp(l)
    GD_error_dd{l} = norm_distance(y0, 0*y0);
    for iter = 1:10
        tSignal = importdata(['./results/forward/FB/forwardSignal_GD_tau', GD.tau{l}, '_lambda', GD.lambda{1}, '_iter', int2str(iter), '.dat'], ' ', 0);
        yi = tSignal(2:end, 1:end-1);
        GD_error_dd{l} = [GD_error_dd{l} norm_distance(y0, yi)];
    end
end

% Plot
x_axis = [0 1 2 3 4 5 6 7 8 9 10];
figure();
colors = winter(length(GD.tau));
for l = 1:length(GD.tau)
    hold on;
    semilogy(x_axis, GD_error_dd{l}, 'Color', colors(l, :), 'Linewidth', 1.5);
end
legend('GD tau 1e16', 'GD tau 2e16', 'GD tau 5e16', 'GD tau 1e17', 'GD tau 2e17', 'GD tau 5e17', 'GD tau 1e18', 'GD tau 2e18');
title(['Dual Distance Error', GD.lambda{1}]);
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;

%====================
% Gradient descent - lambda variation
%====================
disp('GD');
GD.tau    = {'1e18'};
GD.lambda = {'2e-3', '5e-3', '1e-2'};
clear GD_error_dd;
for ii = 1:length(GD.lambda)
    disp(ii)
    GD_error_dd{ii} = norm_distance(y0, 0*y0);
    for iter = 1:10
        tSignal = importdata(['./results/forward/FB/forwardSignal_GD_tau', GD.tau{1}, '_lambda', GD.lambda{ii}, '_iter', int2str(iter), '.dat'], ' ', 0);
        yi = tSignal(2:end, 1:size(y0, 2));
        GD_error_dd{ii} = [GD_error_dd{ii} norm_distance(y0, yi)];
        GD_error_dd{ii}
    end
end

% Plot
x_axis = [0 1 2 3 4 5 6 7 8 9 10];
figure();
axis([0 10 3.6e-17 5e-17])
colors = winter(length(GD.lambda));
for ii = 1:length(GD.lambda)
    hold on;
    semilogy(x_axis, GD_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5);
end
legend('GD lambda 2e-3', 'GD lambda 5e-3', 'GD lambda 1e-2');
title(['Dual Distance Error', GD.tau{1}]);
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
%saveas(gcf, './figures/RD04_GD_dd.fig');

%====================
% Stochastig Gradient descent
%====================
disp('S-GD');
SGD.batch = '3800';
SGD.tau    =  {'6e18'};
SGD.lambda = {'1e-2', '2e-2', '5e-2'};
%SGD.lambda = {'2e-3', '5e-3', '1e-2'};
%SGD.lambda = {'1e-3', '2e-3', '5e-3', '1e-2'};
%SGD.lambda = {'2e-4', '5e-4', '1e-3', '2e-3', '5e-3'};
%SGD.lambda = {'1e-3', '2e-3', '5e-3', '1e-2', '2e-2'};
%SGD.lambda = {'1e-5', '1e-4', '1e-3', '1e-2', '1e-1'};
%SGD.lambda = {'5e-5', '1e-4', '2e-4', '5e-4', '1e-3'};
%SGD.lambda = {'1', '2', '5', '10'};
%SGD.lambda = {'5e-6', '1e-5', '2e-5', '5e-5'};
%SGD.lambda = {'1e-4', '1e-5', '1e-6'};
%SGD.lambda = {'1e-4', '2e-4', '4e-4', '8e-4'};
nIter = 10;
clear SGD_error_dd;
for ii = 1:length(SGD.lambda)
    disp(ii)
    SGD_error_dd{ii} = norm_distance(y0, 0*y0);
    for iter = 1:nIter
        tSignal = importdata(['./results/forward/S-FB/', SGD.batch, '/forwardSignal_S-GD_tau', SGD.tau{1}, '_lambda', SGD.lambda{ii}, '_batch', SGD.batch, '_epoch', int2str(iter), '.dat'], ' ', 0);
        yi = tSignal(2:end, 1:size(y0, 2));
        SGD_error_dd{ii} = [SGD_error_dd{ii} norm_distance(y0, yi)];
        SGD_error_dd{ii}
    end
end
% Plot
x_axis = 0:nIter;
figure();
colors = winter(length(SGD.lambda));
colors_summer = autumn(length(SGD.lambda));
for ii = 1:length(SGD.lambda)
    hold on;
    semilogy(x_axis, SGD_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5);
    %semilogy(0:10, GD_error_dd{ii}, 'Color', colors_summer(ii, :), 'Linewidth', 1.5);
end
axis([0 10 3.6e-17 5e-17])
legend('S-GD 1e-2', 'S-GD 2e-2', 'S-GD 5e-2');
%legend('S-GD 2e-3', 'S-GD 5e-3', 'S-GD 1e-2');
%legend('S-GD 2e-4', 'S-GD 5e-4', 'S-GD 1e-3', 'S-GD 2e-3', 'S-GD 5e-3');
%legend('S-GD lambda 5e-6', 'S-GD lambda 1e-5', 'S-GD lambda 2e-5', 'S-GD lambda 5e-5');
%legend('S-GD lambda 1e-3', 'S-GD lambda 2e-3', 'S-GD lambda 5e-3', 'S-GD lambda 1e-2', 'S-GD lambda 2e-2');
%legend('S-GD tau 1e17', 'S-GD tau 2e17', 'S-GD lambda 5e17');
title(['Dual Distance Error', SGD.tau]);
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;

%====================
% Stochastig Gradient descent - tau variation
%====================
disp('S-GD');
SGD.batch = '2100';
SGD.tau    =  {'2e17', '5e17', '1e18', '2e18'};
%SGD.lambda = {'5e-5', '1e-4', '2e-4', '5e-4', '1e-3'};
SGD.lambda = {'1e-4'};
clear SGD_error_dd;
for ii = 1:length(SGD.tau)
    disp(ii)
    SGD_error_dd{ii} = norm_distance(y0, 0*y0);
    for iter = 1:10
        tSignal = importdata(['./results/forward/S-FB/forwardSignal_S-GD_tau', SGD.tau{ii}, '_lambda', SGD.lambda{1}, '_batch', SGD.batch, '_epoch', int2str(iter), '.dat'], ' ', 0);
        yi = tSignal(2:end, :);
        SGD_error_dd{ii} = [SGD_error_dd{ii} norm_distance(y0, yi)];
        SGD_error_dd{ii}
    end
end
% Plot
x_axis = [0 1 2 3 4 5 6 7 8 9 10];
figure();
colors = winter(length(SGD.tau));
for ii = 1:length(SGD.tau)
    hold on;
    semilogy(x_axis, SGD_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5);
end
axis([0 10 1.15 1.45])
legend('S-GD tau 2e17', 'S-GD tau 5e17', 'S-GD tau 1e18', 'S-GD tau 2e18');
title(['Dual Distance Error - SGD batch 800']);
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
saveas(gcf, ['./figures/RD04_SGD_dd_error_batch800.fig']);



%====================
% FISTA - tau variation
%====================
disp('FISTA');
FISTA.tau    = {'5e16', '1e17', '2e17', '5e17'};
FISTA.lambda = {'1e-2'};
clear FISTA_error_dd;
for l = 1:length(FISTA.tau)
    disp(l)
    FISTA_error_dd{l} = norm_distance(y0, 0*y0);
    for iter = 1:10
        tSignal = importdata(['./results/forward/AFB/forwardSignal_FISTA_tau', FISTA.tau{l}, '_lambda', FISTA.lambda{1}, '_iter', int2str(iter), '.dat'], ' ', 0);
        yi = tSignal(2:end, :);
        FISTA_error_dd{l} = [FISTA_error_dd{l} norm_distance(y0, yi)];
    end
end

% Plot
x_axis = [0 1 2 3 4 5 6 7 8 9 10];
figure();
colors = winter(length(FISTA.tau));
for l = 1:length(FISTA.tau)
    hold on;
    semilogy(x_axis, FISTA_error_dd{l}, 'Color', colors(l, :), 'Linewidth', 1.5);
end
legend('FISTA tau 5e16', 'FISTA tau 1e17', 'FISTA tau 2e17', 'FISTA tau 5e17', 'FISTA tau 1e18', 'FISTA tau 2e18');
title(['Dual Distance Error', FISTA.lambda{1}]);
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;

%====================
% FISTA - lambda variation
%====================
disp('FISTA');
FISTA.tau    = {'2e17'};
FISTA.lambda = {'2e-3', '5e-3', '1e-2'};
clear FISTA_error_dd;
for ii = 1:length(FISTA.lambda)
    disp(ii)
    FISTA_error_dd{ii} = norm_distance(y0, 0*y0);
    for iter = 1:10
        tSignal = importdata(['./results/forward/AFB/forwardSignal_FISTA_tau', FISTA.tau{1}, '_lambda', FISTA.lambda{ii}, '_iter', int2str(iter), '.dat'], ' ', 0);
        yi = tSignal(2:end, 1:size(y0, 2));
        FISTA_error_dd{ii} = [FISTA_error_dd{ii} norm_distance(y0, yi)];
    end
end

% Plot
x_axis = [0 1 2 3 4 5 6 7 8 9 10];
figure();
colors = winter(length(FISTA.lambda));
for ii = 1:length(FISTA.lambda)
    hold on;
    semilogy(x_axis, FISTA_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5);
end
legend('FISTA tau 5e16', 'FISTA tau 1e17', 'FISTA tau 2e17', 'FISTA tau 5e17', 'FISTA tau 1e18', 'FISTA tau 2e18');
title(['Dual Distance Error', FISTA.lambda{ii}]);
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;



%====================
% PDHG - tau variation
%====================
disp('PDHG');
PDHG.sigma = '2';
PDHG.tau    = {'2e17', '5e17', '1e18', '2e18'};
PDHG.lambda = {'1e-3'};
clear PDHG_error_dd;
for ii = 1:length(PDHG.tau)
    disp(ii)
    PDHG_error_dd{ii} = norm_distance(y0, 0*y0);
    for iter = 1:10
        tSignal = importdata(['./results/forward/PDHG/forwardSignal_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau{ii}, '_theta', PDHG.theta, '_lambda', PDHG.lambda{1}, '_iter', int2str(iter), '.dat'], ' ', 0);       
        yi = tSignal(2:end, :);
        PDHG_error_dd{ii} = [PDHG_error_dd{ii} norm_distance(y0, yi)];
    end
end

% Plot
x_axis = [0 1 2 3 4 5 6 7 8 9 10];
figure();
colors = winter(length(PDHG.tau));
for l = 1:2%length(PDHG.tau)
    hold on;
    semilogy(x_axis, PDHG_error_dd{l}, 'Color', colors(l, :), 'Linewidth', 1.5);
end
legend('PDHG tau 2e17', 'PDHG tau 5e17', 'PDHG tau 1e18', 'PDHG tau 2e18');
title(['Dual Distance Error', PDHG.lambda{1}]);
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;


%====================
% PDHG - sigma variation
%====================
disp('PDHG');
PDHG.sigma = {'5e-1', '1'};
PDHG.tau    = {'5e17'};
PDHG.lambda = {'2e-3'};
clear PDHG_error_dd;
for ii = 1:length(PDHG.sigma)
    disp(ii)
    PDHG_error_dd{ii} = norm_distance(y0, 0*y0);
    for iter = 1:10
        tSignal = importdata(['./results/forward/PDHG/forwardSignal_PDHG_sigma', PDHG.sigma{ii}, '_tau', PDHG.tau{1}, '_theta', PDHG.theta, '_lambda', PDHG.lambda{1}, '_iter', int2str(iter), '.dat'], ' ', 0);       
        yi = tSignal(2:end, 1:size(y0, 2));
        PDHG_error_dd{ii} = [PDHG_error_dd{ii} norm_distance(y0, yi)];
    end
end

% Plot
x_axis = [0 1 2 3 4 5 6 7 8 9 10];
figure();
colors = winter(length(PDHG.sigma));
for ii = 1:length(PDHG.sigma)
    hold on;
    plot(x_axis, PDHG_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5);
end
legend('PDHG sigma 5e-1', 'PDHG sigma 1');
title(['Dual Distance Error lambda', PDHG.lambda{1}]);
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;



%====================
% SPDHG - lambda variation
%====================
disp('SPDHG');
SPDHG.sigma = {'1'};
SPDHG.tau    = {'2e18'};
SPDHG.lambda = {'2e-3', '5e-3', '1e-2'};
%SPDHG.lambda = {'5e-4', '1e-3'};
%SPDHG.lambda = {'5e-4', '2e-3', '5e-3', '1e-2'};
SPDHG.batch = '2100';
clear SPDHG_error_dd;
for ii = 1:length(SPDHG.lambda)
    disp(ii)
    SPDHG_error_dd{ii} = norm_distance(y0, 0*y0);
    for iter = 1:9
        tSignal = importdata(['./results/forward/S-PDHG/forwardSignal_S-PDHG_sigma', SPDHG.sigma{1}, '_tau', SPDHG.tau{1}, '_theta', SPDHG.theta, '_lambda', SPDHG.lambda{ii}, '_batch', SPDHG.batch, '_epoch', int2str(iter), '.dat'], ' ', 0);       
        yi = tSignal(2:end, :);
        SPDHG_error_dd{ii} = [SPDHG_error_dd{ii} norm_distance(y0, yi)];
        SPDHG_error_dd{ii}
    end
end

% Plot
x_axis = [0 1 2 3 4 5 6 7 8 9];
figure();
colors = winter(length(SPDHG.lambda));
for ii = 1:length(SPDHG.lambda)
    hold on;
    plot(x_axis, SPDHG_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5);
end
%axis([0 10 1.15 1.45])
legend('2e-3', '5e-3', '1e-2');
title(['Dual Distance Error tau ', SPDHG.tau{1}]);
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;

%************************************************************************************************************************
%========================================================================================================================
% PRIMAL ERROR
%========================================================================================================================
%************************************************************************************************************************
% Load Initial Pressure
u0Matrix = importdata('./input_data/reconstruction_15iterations.dat', ' ', 0);
u0 = matrix2cube(u0Matrix, Nz);
h = plot_projection(u0, dx);


%====================
% Gradient descent - lambda variation
%====================
disp('GD');
GD.tau    = {'1e18'};
GD.lambda = {'2e-3', '5e-3', '1e-2'};
clear GD_error_pd;
for ii = 1:length(GD.lambda)
    disp(ii)
    GD_error_pd{ii} = norm_distance(u0, 0*u0);
    for iter = 1:10
        pixelPressure = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau{1}, '_lambda', GD.lambda{ii}, '_iter', int2str(iter), '.dat'], ' ', 0);
        pixelPressure = max(0, matrix2cube(pixelPressure, Nz));
        GD_error_pd{ii} = [GD_error_pd{ii} norm_distance(u0, pixelPressure)];
    end
end

% Plot
x_axis = [0 1 2 3 4 5 6 7 8 9 10];
figure();
%axis([0 10 3.6e-17 5e-17])
colors = winter(length(GD.lambda));
for ii = 1:length(GD.lambda)
    semilogy(x_axis, GD_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5);
    hold on;
end
axis([0 10 1 100]);
legend('GD lambda 2e-3', 'GD lambda 5e-3', 'GD lambda 1e-2');
title(['Primal Distance Error', GD.tau{1}]);
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
%saveas(gcf, './figures/RD04_GD_dd.fig');

%====================
% Stochastic Gradient descent - lambda variation
%====================
disp('S-GD');
SGD.batch = '3800';
SGD.tau    =  {'2e18'};
SGD.lambda = {'1e-3', '2e-3', '5e-3', '1e-2', '2e-2'};
clear SGD_error_pd;
for ii = 1:length(SGD.lambda)
    disp(ii)
    SGD_error_pd{ii} = norm_distance(u0, 0*u0);
    for iter = 1:10
        pixelPressure = importdata(['./results/adjoint/S-FB/pixelPressure_S-GD_tau', SGD.tau{1}, '_lambda', SGD.lambda{ii}, '_batch', SGD.batch, '_epoch', int2str(iter), '.dat'], ' ', 0);
        pixelPressure = max(0, matrix2cube(pixelPressure, Nz));
        SGD_error_pd{ii} = [SGD_error_pd{ii} norm_distance(u0, pixelPressure)];
        SGD_error_pd{ii}
    end
end
% Plot
x_axis = [0 1 2 3 4 5 6 7 8 9 10];
figure();
colors = winter(length(SGD.lambda));
for ii = 1:length(SGD.lambda)
    semilogy(x_axis, SGD_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5);
    hold on;
end
axis([0 10 1 200]);
legend('S-GD 1e-3', 'S-GD 2e-3', 'S-GD 5e-3', 'S-GD 1e-2', 'S-GD 2e-2');
title(['Dual Distance Error', SGD.tau]);
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;



%========================================================================================================================
% LONG RECONSTRUCTION
%========================================================================================================================
%====================
% Gradient descent - tau variation
%====================
nIter = 20;
disp('GD');
GD.tau    = '5e17';
GD.lambda = '5e-3';
clear GD_error_dd;
GD_error_dd = norm_distance(y0, 0*y0);
for iter = 1:nIter
    tSignal = importdata(['./results/forward/FB/forwardSignal_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
    yi = tSignal(2:end, 1:end);
    GD_error_dd = [GD_error_dd norm_distance(y0, yi)];
    GD_error_dd(end)
end

% Plot
x_axis = 0:nIter;
figure();
hold on;
semilogy(x_axis, GD_error_dd, 'Linewidth', 1.5);
%legend('GD tau 1e16', 'GD tau 2e16', 'GD tau 5e16', 'GD tau 1e17', 'GD tau 2e17', 'GD tau 5e17', 'GD tau 1e18', 'GD tau 2e18');
title(['Dual Distance Error']);
grid on;
box on;
axis([0 20 3.5e-17 5e-17])
ax = gca;
ax.GridAlpha = 0.2;


%==============================
% Gradient Descent
%==============================
GD.iter = '15';
pixelPressureMatrix = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', GD.iter, '.dat'], ' ', 0);
pixelPressure_recon = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure_recon, dx);
a = axes;
t = title(['GD - tau = ', GD.tau, ', lambda = ', GD.lambda, ', iter = ', GD.iter]);
a.Visible = 'off'; 
t.Visible = 'on';
recon_15iter = cube2matrix(pixelPressure_recon);
dlmwrite('./input_data/reconstruction_15iterations.dat', recon_15iter, ' ');
%saveas(gcf, ['./figures/RD06_GD_tau', GD.tau{i}, '_lambda', GD.lambda{j}, '_iter', GD.iter, '.fig']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('GD');
GD.tau    = {'5e17'};
GD.lambda = {'5e-3'};
clear GD_error_pd;
for ii = 1:length(GD.lambda)
    disp(ii)
    GD_error_pd{ii} = norm_distance(u0, 0*u0);
    for iter = 1:15
        pixelPressure = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau{1}, '_lambda', GD.lambda{ii}, '_iter', int2str(iter), '.dat'], ' ', 0);
        pixelPressure = max(0, matrix2cube(pixelPressure, Nz));
        GD_error_pd{ii} = [GD_error_pd{ii} norm_distance(u0, pixelPressure)];
    end
end

% Plot
x_axis = 0:15;
figure();
%axis([0 10 3.6e-17 5e-17])
colors = winter(length(GD.lambda));
for ii = 1:length(GD.lambda)
    semilogy(x_axis, GD_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5);
    hold on;
end
title(['Primal Distance Error', GD.tau{1}]);
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
%saveas(gcf, './figures/RD04_GD_dd.fig');

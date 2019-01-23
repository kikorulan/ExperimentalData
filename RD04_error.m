% Read data from files
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/RD04_palm;

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

%========================================================================================================================
% DUAL DATA
%========================================================================================================================
% Import data
filenameData = './input_data/forwardSignal_reference_8736sensors.dat';
timeSignal = importdata(filenameData, ' ', 0);
timeRD = timeSignal(1, :);
y0 = timeSignal(2:end, :);
% Plot
figure;
imagesc(y0);
box on;
colorbar();

%==================================================
% Adjoint
%==================================================
% RT
filenameData = './input_data/pixelPressure_adjoint_8736sensors.dat';
pixelPressureMatrix = importdata(filenameData, ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
saveas(gcf, './figures/RTadjoint', 'png');

% KWave
filenameData = './input_data/pressure_adjoint_kWave_8736sensors.dat';
pixelPressureMatrix = importdata(filenameData, ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
saveas(gcf, './figures/KWadjoint', 'png');


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
pixelPressureMatrix = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', GD.iter, '.dat'], ' ', 0);
pixelPressure_GD = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure_GD, dx);
a = axes;
a.Visible = 'off'; 
t = title(['GD - t = ', GD.tau, ', l = ', GD.lambda, ', iter = ', GD.iter, ' - homogeneous SS']);
t.Visible = 'on'; 
saveas(gcf, ['./figures/RD04_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', GD.iter, '.fig']);
saveas(gcf, ['./figures/RD04_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', GD.iter], 'epsc');

%==============================
% Stochastic Gradient Descent
%==============================
pixelPressureMatrix = importdata(['./results/adjoint/S-FB/pixelPressure_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_epoch', SGD.epoch, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['S-GD - t = ', SGD.tau, ', l = ', SGD.lambda, ', batch = ', SGD.batch, ', epoch = ', SGD.epoch, ' - homogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, ['./figures/RD04_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_epoch', SGD.epoch, '.fig']);

%==============================
% FISTA
%==============================
pixelPressureMatrix = importdata(['./results/adjoint/AFB/pixelPressure_FISTA_tau', FISTA.tau, '_lambda', FISTA.lambda, '_iter', FISTA.iter, '.dat'], ' ', 0);
pixelPressure_FISTA = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure_FISTA, dx);
a = axes;
a.Visible = 'off'; 
t = title(['FISTA - t = ', FISTA.tau, ', l = ', FISTA.lambda, ', iter = ', FISTA.iter, ' - homogeneous SS']);
t.Visible = 'on'; 
saveas(gcf, ['./figures/RD04_FISTA_tau', FISTA.tau, '_lambda', FISTA.lambda, '_iter', FISTA.iter, '.fig']);
saveas(gcf, ['./figures/RD04_FISTA_tau', FISTA.tau, '_lambda', FISTA.lambda, '_iter', FISTA.iter], 'epsc');

%==============================
% PDHG
%==============================
pixelPressureMatrix = importdata(['./results/adjoint/PDHG/pixelPressure_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau, '_theta', PDHG.theta, '_lambda', PDHG.lambda, '_iter', PDHG.iter, '.dat'], ' ', 0);
pixelPressure_PDHG = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure_PDHG, dx);
a = axes;
a.Visible = 'off'; 
t = title(['PDHG - s = ', PDHG.sigma, ', t = ', PDHG.tau, ', l = ', PDHG.lambda, ', iter = ', PDHG.iter, ' - homogeneous SS']);
t.Visible = 'on'; 
saveas(gcf, ['./figures/RD04_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau, '_lambda', PDHG.lambda, '_iter', PDHG.iter, '.fig']);
saveas(gcf, ['./figures/RD04_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau, '_lambda', PDHG.lambda, '_iter', PDHG.iter], 'epsc');

%==============================
% S-PDHG
%==============================
pixelPressureMatrix = importdata(['./results/adjoint/S-PDHG/pixelPressure_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau, '_theta', SPDHG.theta, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_epoch', SPDHG.epoch, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['S-PDHG - s = ', SPDHG.sigma, ', t = ', SPDHG.tau, ', l = ', SPDHG.lambda, ', batch = ', SPDHG.batch, ', epoch = ', SPDHG.epoch, ' - homogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, ['./figures/RD04_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_epoch', SPDHG.epoch, '.fig']);

%========================================================================================================================
% DUAL DISTANCE ERROR
%========================================================================================================================
disp('******* DUAL DISTANCE ********');
% Gradient descent
disp('GD');
GD_error_dd = norm_distance(y0, 0*y0);
for iter = 1:10
    tSignal = importdata(['./results/forward/FB/forwardSignal_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
    yi = tSignal(2:end, :);
    GD_error_dd = [GD_error_dd norm_distance(y0, yi)];
end

% Stochastic Gradient descent
disp('S-GD');
SGD_error_dd = norm_distance(y0, 0*y0);
for iter = 1:10
    tSignal = importdata(['./results/forward/S-FB/forwardSignal_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_epoch', int2str(iter), '.dat'], ' ', 0);
    yi = tSignal(2:end, :);
    SGD_error_dd = [SGD_error_dd norm_distance(y0, yi)];
end

% FISTA
disp('FISTA');
FISTA_error_dd = norm_distance(y0, 0*y0);
for iter = 1:10
    tSignal = importdata(['./results/forward/AFB/forwardSignal_FISTA_tau', FISTA.tau, '_lambda', FISTA.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
    yi = tSignal(2:end, :);
    FISTA_error_dd = [FISTA_error_dd norm_distance(y0, yi)];
end

% PDHG
disp('PDHG');
PDHG_error_dd = norm_distance(y0, 0*y0);
for iter = 1:10
    tSignal = importdata(['./results/forward/PDHG/forwardSignal_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau, '_theta', PDHG.theta, '_lambda', PDHG.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
    yi = tSignal(2:end, :);
    PDHG_error_dd = [PDHG_error_dd norm_distance(y0, yi)];
end

% SPDHG
disp('S-PDHG');
SPDHG.tau    = '5e17';
SPDHG_error_dd = norm_distance(y0, 0*y0);
for iter = 1:10
    tSignal = importdata(['./results/forward/S-PDHG/forwardSignal_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau, '_theta', SPDHG.theta, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_epoch', int2str(iter), '.dat'], ' ', 0); 
    yi = tSignal(2:end, :);
    SPDHG_error_dd = [SPDHG_error_dd norm_distance(y0, yi)];
end

% Plot
x_axis = [0 1 2 3 4 5 6 7 8 9 10];
figure();
plot(x_axis, GD_error_dd, 'Color', 'r', 'Linewidth', 1.5);
hold on;
plot(x_axis, SGD_error_dd, 'Color', 'g', 'Linewidth', 1.5);
plot(x_axis, FISTA_error_dd, 'Color', 'b', 'Linewidth', 1.5);
plot(x_axis, PDHG_error_dd, 'Color', 'm', 'Linewidth', 1.5);
plot(x_axis, SPDHG_error_dd, 'Color', 'c', 'Linewidth', 1.5);
legend('GD', 'S-GD', 'FISTA', 'PDHG', 'S-PDHG');
%legend('GD', 'FISTA', 'PDHG', 'S-PDHG');
%legend('FB', 'AFB', 'PDHG');
title('Dual Distance Error');
axis([0 10 1.15 1.45])
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
saveas(gcf, ['./figures/RD04_dd_error.fig']);
saveas(gcf, ['./figures/RD04_dd_error'], 'epsc');

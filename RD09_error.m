% Read data from files
%cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/RD09_finger2_doubleRes_subsampled;
cd /scratch0/NOT_BACKED_UP/frullan//ExperimentalData/RD09_finger2_doubleRes_subsampled;

clear all;
close all;

% Functions
norm_distance = @(x, y) sqrt(sum((x(:) - y(:)).*(x(:) - y(:))));

%==================================================
% Dimensions
%==================================================
% Import dimensions
dim = importdata('./input_data/dimensions.dat', ' ', 0);
Nx = dim(1, 1); dx = dim(2, 1);
Ny = dim(1, 2); dy = dim(2, 2);
Nz = dim(1, 3); dz = dim(2, 3);

%========================================================================================================================
% PRIMAL AND DUAL DATA
%========================================================================================================================
% Forward signal
time_signal = importdata(['./input_data/forwardSignal_reference_3600sensors.dat'], ' ', 0);
y0 = time_signal(2:end, :);

% Load Adjoint Pressure
pressure_adjoint = importdata('./input_data/pixelPressure_adjoint_3600sensors.dat', ' ', 0);
pressure_adjoint = matrix2cube(pressure_adjoint, Nz);
h = plot_projection(pressure_adjoint, dx);
a = axes;
t = title('Adjoint RT');
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, './figures/RD09_adjointRT.fig');

%========================================================================================================================
% ITERATIVE RECONSTRUCTION
%========================================================================================================================
iter = 5;

%==============================
% Gradient Descent
%==============================
% GD **************************
GD = [];
GD.tau    = '2e18';
GD.lambda = '2e-3';
GD.iter   = int2str(iter);
%******************************
pixelPressureMatrix = importdata(['./results_without_zero_layer/adjoint/FB/pixelPressure_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', GD.iter, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['FB - t = ', GD.tau, ', l = ', GD.lambda, ', iter = ', GD.iter, ' - homogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, ['./figures/RD09_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', GD.iter, '.fig']);

%==============================
% Stochastic Gradient Descent
%==============================
% S-GD ************************
SGD = [];
SGD.tau    = '4e18';
SGD.lambda = '1e-3';
SGD.batch  = '600';
SGD.epoch  = int2str(iter);
%******************************
pixelPressureMatrix = importdata(['./results_without_zero_layer/adjoint/S-FB/', SGD.batch, '/pixelPressure_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_subepoch', SGD.epoch, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['S-FB - t = ', SGD.tau, ', l = ', SGD.lambda, ', batch = ', SGD.batch, ', epoch = ', SGD.epoch, ' - homogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, ['./figures/RD09_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_epoch', SGD.epoch, '.fig']);

%==============================
% FISTA
%==============================
% FISTA ***********************
FISTA = [];
FISTA.tau    = '1e18';
FISTA.lambda = '5e-3';
FISTA.iter   = int2str(iter);
%******************************
pixelPressureMatrix = importdata(['./results_without_zero_layer/adjoint/AFB/pixelPressure_FISTA_tau', FISTA.tau, '_lambda', FISTA.lambda, '_iter', FISTA.iter, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['AFB - t = ', FISTA.tau, ', l = ', FISTA.lambda, ', iter = ', FISTA.iter, ' - homogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, ['./figures/RD09_FISTA_tau', FISTA.tau, '_lambda', FISTA.lambda, '_iter', FISTA.iter, '.fig']);

%==============================
% PDHG
%==============================
% PDHG ************************
PDHG = [];
PDHG.sigma  = '1';
PDHG.tau    = '2e18';
PDHG.theta  = '1';
PDHG.lambda = '2e-3';
PDHG.iter   = int2str(iter);
%******************************
pixelPressureMatrix = importdata(['./results_without_zero_layer/adjoint/PDHG/pixelPressure_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau, '_theta', PDHG.theta, '_lambda', PDHG.lambda, '_iter', PDHG.iter, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['PDHG - s = ', PDHG.sigma, ', t = ', PDHG.tau, ', l = ', PDHG.lambda, ', iter = ', PDHG.iter, ' - homogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, ['./figures/RD09_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau, '_lambda', PDHG.lambda, '_iter', PDHG.iter, '.fig']);

%==============================
% S-PDHG
%==============================
% S-PDHG **********************
SPDHG = [];
SPDHG.sigma  = '1';
SPDHG.tau    = '2e18';
SPDHG.theta  = '1';
SPDHG.lambda = '2e-3';
SPDHG.batch  = '600';
SPDHG.epoch  = int2str(iter);
%******************************
pixelPressureMatrix = importdata(['./results_without_zero_layer/adjoint/S-PDHG/', SPDHG.batch, '/pixelPressure_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau, '_theta', SPDHG.theta, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_subepoch', SPDHG.epoch, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['S-PDHG - s = ', SPDHG.sigma, ', t = ', SPDHG.tau, ', l = ', SPDHG.lambda, ', batch = ', SPDHG.batch, ', epoch = ', SPDHG.epoch, ' - homogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, ['./figures/RD09_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_epoch', SPDHG.epoch, '.fig']);


%========================================================================================================================
% DUAL DISTANCE ERROR
%========================================================================================================================
disp('******* DUAL DISTANCE ********');
nIter = 5;
x_axis = 0:nIter;
NORM = 2.5e17;

% Gradient descent
disp('GD');
GD.lambda = '2e-3';
GD.tau = '2e18';
clear GD_error_dd;
GD_error_dd = norm_distance(y0, 0*y0);
for iter = 1:5
    tSignal = importdata(['./results_without_zero_layer/forward/FB/forwardSignal_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
    yi = tSignal(2:end, :);
    GD_error_dd = [GD_error_dd norm_distance(y0, NORM*yi)];
end

% Stochastic Gradient descent
disp('S-GD');
SGD.lambda = '1e-3';
SGD.batch = '600';
SGD.tau = '4e18';
clear SGD_error_dd;
SGD_error_dd = norm_distance(y0, 0*y0);
for iter = 1:5
    tSignal = importdata(['./results_without_zero_layer/forward/S-FB/', SGD.batch, '/forwardSignal_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_subepoch', int2str(iter), '.dat'], ' ', 0);
    yi = tSignal(2:end, :);
    SGD_error_dd = [SGD_error_dd norm_distance(y0, NORM*yi)];
end


% FISTA
disp('FISTA');
FISTA.lambda = '5e-3';
FISTA.tau =  '1e18';
clear FISTA_error_dd;
FISTA_error_dd = norm_distance(y0, 0*y0);
for iter = 1:5
    tSignal = importdata(['./results_without_zero_layer/forward/AFB/forwardSignal_FISTA_tau', FISTA.tau, '_lambda', FISTA.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
    yi = tSignal(2:end, :);
    FISTA_error_dd = [FISTA_error_dd norm_distance(y0, NORM*yi)];
end

% PDHG
disp('PDHG');
PDHG.lambda = '2e-3';
PDHG.tau = '2e18';
clear PDHG_error_dd;
PDHG_error_dd = norm_distance(y0, 0*y0);
for iter = 1:5
    tSignal = importdata(['./results_without_zero_layer/forward/PDHG/forwardSignal_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau, '_theta', PDHG.theta, '_lambda', PDHG.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
    yi = tSignal(2:end, :);
    PDHG_error_dd = [PDHG_error_dd norm_distance(y0, NORM*yi)];
end


% SPDHG
disp('S-PDHG');
SPDHG.lambda = '2e-3';
SPDHG.tau = '4e18';
SPDHG.batch = '600';
clear SPDHG_error_dd;
SPDHG_error_dd = norm_distance(y0, 0*y0);
for iter = 1:5
    tSignal = importdata(['./results_without_zero_layer/forward/S-PDHG/', SPDHG.batch, '/forwardSignal_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau, '_theta', SPDHG.theta, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_subepoch', int2str(iter), '.dat'], ' ', 0); 
    yi = tSignal(2:end, :);
    SPDHG_error_dd = [SPDHG_error_dd norm_distance(y0, NORM*yi)];
end

%==============================
% PLOT ALL
%==============================
% Plot
figure();
semilogy(x_axis, GD_error_dd, 'Color', 'r', 'Linewidth', 1.5);
hold on;
semilogy(x_axis, SGD_error_dd, 'Color', 'g', 'Linewidth', 1.5);
semilogy(x_axis, FISTA_error_dd, 'Color', 'b', 'Linewidth', 1.5);
semilogy(x_axis, PDHG_error_dd, 'Color', 'm', 'Linewidth', 1.5);
semilogy(x_axis, SPDHG_error_dd, 'Color', 'c', 'Linewidth', 1.5);
legend('GD', 'S-GD', 'FISTA', 'PDHG', 'S-PDHG');
title('Dual Distance Error - homogeneous SS');
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 nIter 3.5 6])
saveas(gcf, ['./figures/RD09_dd_error.fig']);


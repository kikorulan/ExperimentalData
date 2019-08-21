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
figure;
imagesc(y0);
% Forward signal
time_signal = importdata(['./input_data/forwardSignal_reference_norm_3600sensors.dat'], ' ', 0);
y0_norm = time_signal(2:end, :);

% Load Adjoint Pressure
pressure_adjoint = importdata('./input_data/pixelPressure_adjoint_3600sensors.dat', ' ', 0);
pressure_adjoint = matrix2cube(pressure_adjoint, Nz);
h = plot_projection(pressure_adjoint, dx);

% Display maximum
disp(['Max y0: ', num2str(max(y0(:)))])
disp(['Max y0 norm: ', num2str(max(y0_norm(:)))])
disp(['Max pressure adjoint: ', num2str(max(pressure_adjoint(:)))])

%========================================================================================================================
% ITERATIVE RECONSTRUCTION
%========================================================================================================================
iter = 5;

%==============================
% Gradient Descent
%==============================
% GD **************************
GD = [];
GD.tau    = '4e18';
GD.lambda = '5e-3';
GD.iter   = int2str(iter);
%******************************
pixelPressureMatrix = importdata(['./results_with_zero_layer/adjoint/FB/pixelPressure_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', GD.iter, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['GD - t = ', GD.tau, ', l = ', GD.lambda, ', iter = ', GD.iter, ' - homogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
%saveas(gcf, ['./figures/Example80_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', GD.iter, '.fig']);

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
pixelPressureMatrix = importdata(['./results/adjoint/S-FB/', SGD.batch, '/pixelPressure_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_subepoch', SGD.epoch, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['S-GD - t = ', SGD.tau, ', l = ', SGD.lambda, ', batch = ', SGD.batch, ', epoch = ', SGD.epoch, ' - homogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
%saveas(gcf, ['./figures/Example80_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_epoch', SGD.epoch, '.fig']);

%==============================
% FISTA
%==============================
% FISTA ***********************
FISTA = [];
FISTA.tau    = '2e17';
FISTA.lambda = '1e-3';
FISTA.iter   = int2str(iter);
%******************************
%pixelPressureMatrix = importdata(['./results/adjoint/AFB/pixelPressure_FISTA_tau', FISTA.tau, '_lambda', FISTA.lambda, '_iter', FISTA.iter, '.dat'], ' ', 0);
pixelPressureMatrix = importdata(['./results/ground_truth/pixelPressure_FISTA_tau', FISTA.tau, '_lambda', FISTA.lambda, '_iter', FISTA.iter, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['FISTA - t = ', FISTA.tau, ', l = ', FISTA.lambda, ', iter = ', FISTA.iter, ' - homogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
%saveas(gcf, ['./figures/Example80_FISTA_tau', FISTA.tau, '_lambda', FISTA.lambda, '_iter', FISTA.iter, '.fig']);

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
pixelPressureMatrix = importdata(['./results/adjoint/PDHG/pixelPressure_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau, '_theta', PDHG.theta, '_lambda', PDHG.lambda, '_iter', PDHG.iter, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['PDHG - s = ', PDHG.sigma, ', t = ', PDHG.tau, ', l = ', PDHG.lambda, ', iter = ', PDHG.iter, ' - homogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
%saveas(gcf, ['./figures/Example80_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau, '_lambda', PDHG.lambda, '_iter', PDHG.iter, '.fig']);

%==============================
% S-PDHG
%==============================
% S-PDHG **********************
SPDHG = [];
SPDHG.sigma  = '1';
SPDHG.tau    = '2e18';
SPDHG.theta  = '1';
SPDHG.lambda = '1e-3';
SPDHG.batch  = '200';
SPDHG.epoch  = int2str(iter);
%******************************
pixelPressureMatrix = importdata(['./results/adjoint/S-PDHG/', SPDHG.batch, '/pixelPressure_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau, '_theta', SPDHG.theta, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_subepoch', SPDHG.epoch, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['S-PDHG - s = ', SPDHG.sigma, ', t = ', SPDHG.tau, ', l = ', SPDHG.lambda, ', batch = ', SPDHG.batch, ', epoch = ', SPDHG.epoch, ' - homogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
%saveas(gcf, ['./figures/Example80_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_epoch', SPDHG.epoch, '.fig']);


%========================================================================================================================
% DUAL DISTANCE ERROR
%========================================================================================================================
disp('******* DUAL DISTANCE ********');
nIter = 5;
x_axis = 0:nIter;
NORM = 2.5e17;
%==============================
% Gradient descent
%==============================
disp('GD');
GD.lambda = '2e-3';
%GD.tau = {'1e18', '2e18'};
GD.tau = {'1e18', '2e18', '4e18'};
L = length(GD.tau);
clear GD_error_dd;
for ii = 1:L
    disp(ii)
    GD_error_dd{ii} = norm_distance(y0, 0*y0);
    for iter = 1:5
        tSignal = importdata(['./results/forward/FB/forwardSignal_GD_tau', GD.tau{ii}, '_lambda', GD.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
        yi = tSignal(2:end, :);
        GD_error_dd{ii} = [GD_error_dd{ii} norm_distance(y0, NORM*yi)];
    end
end
% Plot
figure();
colors = winter(L);
for ii = 1:L
    semilogy(x_axis, GD_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 nIter 3 6]);
%legend('GD 1e18', 'GD 2e18');
legend('GD 1e18', 'GD 2e18', 'GD 4e18');

%==============================
% Stochastic Gradient descent
%==============================
disp('S-GD');
SGD.lambda = '1e-3';
SGD.batch = '600'
SGD.tau = {'1e18', '2e18', '4e18'};
L = length(SGD.tau);
clear SGD_error_dd;
for ii = 1:L
    disp(ii)
    SGD_error_dd{ii} = norm_distance(y0, 0*y0);
    for iter = 1:5
        tSignal = importdata(['./results/forward/S-FB/', SGD.batch, '/forwardSignal_S-GD_tau', SGD.tau{ii}, '_lambda', SGD.lambda, '_batch', SGD.batch, '_subepoch', int2str(iter), '.dat'], ' ', 0);
        yi = tSignal(2:end, :);
        SGD_error_dd{ii} = [SGD_error_dd{ii} norm_distance(y0, NORM*yi)];
    end
end
% Plot
figure();
colors = winter(L);
for ii = 1:L
    semilogy(x_axis, SGD_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 nIter 3 6]);
legend('SGD 1e18', 'SGD 2e18', 'SGD 4e18');


%==============================
% FISTA
%==============================
disp('FISTA');
FISTA.lambda = '5e-3';
FISTA.tau =  {'1e18', '2e18', '4e18'};
L = length(FISTA.tau);
clear FISTA_error_dd;
for ii = 1:L
    disp(ii)
    FISTA_error_dd{ii} = norm_distance(y0, 0*y0);
    for iter = 1:5
        tSignal = importdata(['./results/forward/AFB/forwardSignal_FISTA_tau', FISTA.tau{ii}, '_lambda', FISTA.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
        yi = tSignal(2:end, :);
        FISTA_error_dd{ii} = [FISTA_error_dd{ii} norm_distance(y0, NORM*yi)];
    end
end
% Plot
figure();
colors = winter(L);
for ii = 1:L
    semilogy(x_axis, FISTA_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 nIter 3 6]);
legend('FISTA 1e18', 'FISTA 2e18', 'FISTA 4e18');

%==============================
% PDHG
%==============================
disp('PDHG');
PDHG.lambda = '2e-3';
PDHG.tau = {'1e18', '2e18', '4e18'};
L = length(PDHG.tau);
clear PDHG_error_dd;
for ii = 1:L
    disp(ii)
    PDHG_error_dd{ii} = norm_distance(y0, 0*y0);
    for iter = 1:5
        tSignal = importdata(['./results/forward/PDHG/forwardSignal_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau{ii}, '_theta', PDHG.theta, '_lambda', PDHG.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
        yi = tSignal(2:end, :);
        PDHG_error_dd{ii} = [PDHG_error_dd{ii} norm_distance(y0, NORM*yi)];
    end
end
% Plot
figure();
colors = winter(L);
for ii = 1:L
    semilogy(x_axis, PDHG_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 nIter 3 6]);
legend('PDHG 1e18', 'PDHG 2e18', 'PDHG 4e18');


%==============================
% SPDHG
%==============================
disp('S-PDHG');
SPDHG.lambda = '5e-3';
SPDHG.tau = {'2e18', '4e18', '8e18'};
SPDHG.batch = '1200';
L = length(SPDHG.tau);
clear SPDHG_error_dd;
for ii = 1:L
    disp(ii)
    SPDHG_error_dd{ii} = norm_distance(y0, 0*y0);
    for iter = 1:5
        tSignal = importdata(['./results/forward/S-PDHG/', SPDHG.batch, '/forwardSignal_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau{ii}, '_theta', SPDHG.theta, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_subepoch', int2str(iter), '.dat'], ' ', 0); 
        yi = tSignal(2:end, :);
        SPDHG_error_dd{ii} = [SPDHG_error_dd{ii} norm_distance(y0, NORM*yi)];
    end
end
% Plot
figure();
colors = winter(L);
for ii = 1:L
    semilogy(x_axis, SPDHG_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 nIter 3 6]);
legend('SPDHG 2e18', 'SPDHG 4e18', 'SPDHG 8e18');


%==============================
% PLOT ALL
%==============================
% Plot
figure();
semilogy(x_axis, GD_error_dd{2}, 'Color', 'r', 'Linewidth', 1.5);
hold on;
semilogy(x_axis, SGD_error_dd{3}, 'Color', 'g', 'Linewidth', 1.5);
semilogy(x_axis, FISTA_error_dd{2}, 'Color', 'b', 'Linewidth', 1.5);
semilogy(x_axis, PDHG_error_dd{2}, 'Color', 'm', 'Linewidth', 1.5);
%semilogy(x_axis, SPDHG_error_dd{1}, 'Color', 'c', 'Linewidth', 1.5);
legend('GD', 'S-GD', 'FISTA', 'PDHG', 'S-PDHG');
title('Dual Distance Error - homogeneous SS');
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 nIter 3 6])
%saveas(gcf, ['./figures/Example80_dd_error.fig']);


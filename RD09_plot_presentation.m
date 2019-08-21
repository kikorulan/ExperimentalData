% Read data from files
cd /scratch0/NOT_BACKED_UP/frullan/ExperimentalData/RD09_finger2_doubleRes_subsampled;

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
imagesc(y0(1:60, :));
h = colorbar();
set(gca,'FontSize', 20);
yt = get(h,'YTick');
exponent = 2;
h.Ruler.Exponent = -2;
for i = 1:length(h.TickLabels)
    h.TickLabels{i} = sprintf('%2.1f',(10^exponent)*yt(i));
end
xlabel('time step')
ylabel('sensor')
%saveas(gcf, './figures_paper/RD09_forwardSignal', 'epsc');

% Load Adjoint Pressure
pressure_adjoint = importdata('./input_data/pixelPressure_adjoint_3600sensors.dat', ' ', 0);
pressure_adjoint = matrix2cube(pressure_adjoint, Nz);
h = plot_projection_compact(pressure_adjoint, dx, false);
%saveas(gcf, './figures_paper/RD09_pressureAdjoint', 'epsc');

%========================================================================================================================
% ITERATIVE RECONSTRUCTION
%========================================================================================================================

%==============================
% Gradient Descent
%==============================
% GD **************************
GD = [];
GD.tau    = '8e17';
GD.lambda = '1e-4';
GD.iter   = int2str(20);
%******************************
pixelPressureMatrix = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', GD.iter, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection_compact(pixelPressure, dx, false);
%saveas(gcf, './figures_paper/RD09_GD', 'epsc');

%==============================
% Stochastic Gradient Descent
%==============================
% S-GD ************************
SGD = [];
SGD.tau    = '1.6e18';
SGD.lambda = '1e-4';
SGD.batch  = '100';
SGD.epoch  = int2str(20);
%******************************
pixelPressureMatrix = importdata(['./results/adjoint/S-FB/pixelPressure_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_subepoch', SGD.epoch, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection_compact(pixelPressure, dx, false);
%saveas(gcf, './figures_paper/RD09_S-GD', 'epsc');

%==============================
% FISTA
%==============================
% FISTA ***********************
FISTA = [];
FISTA.tau    = '2e17';
FISTA.lambda = '1e-4';
FISTA.iter   = int2str(20);
%******************************
pixelPressureMatrix = importdata(['./results/adjoint/AFB/pixelPressure_FISTA_tau', FISTA.tau, '_lambda', FISTA.lambda, '_iter', FISTA.iter, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection_compact(pixelPressure, dx, false);
%saveas(gcf, './figures_paper/RD09_FISTA', 'epsc');

%==============================
% PDHG
%==============================
% PDHG ************************
PDHG = [];
PDHG.sigma  = '1';
PDHG.tau    = '8e17';
PDHG.theta  = '1';
PDHG.lambda = '1e-4';
PDHG.iter   = int2str(20);
%******************************
pixelPressureMatrix = importdata(['./results/adjoint/PDHG/pixelPressure_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau, '_theta', PDHG.theta, '_lambda', PDHG.lambda, '_iter', PDHG.iter, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection_compact(pixelPressure, dx, false);
%saveas(gcf, './figures_paper/RD09_PDHG', 'epsc');


%==============================
% S-PDHG
%==============================
% S-PDHG **********************
SPDHG = [];
SPDHG.sigma  = '1';
SPDHG.tau    = '5e16';
SPDHG.theta  = '1';
SPDHG.lambda = '1e-4';
SPDHG.batch  = '100';
SPDHG.epoch  = int2str(20);
%******************************
pixelPressureMatrix = importdata(['./results/adjoint/S-PDHG/pixelPressure_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau, '_theta', SPDHG.theta, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_subepoch', SPDHG.epoch, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection_compact(pixelPressure, dx, false);
saveas(gcf, './figures_paper/RD09_S-PDHG', 'epsc');

%==========================================================================================================================================================================
%===============================                                                            ===============================================================================
%===============================                     DISTANCE ERROR                         ===============================================================================
%===============================                                                            ===============================================================================
%==========================================================================================================================================================================
disp('******* DUAL DISTANCE ********');
NORM = 2.5e17;

%======================================================================
% PLOT ALL DUAL
%======================================================================
load ./results/error_vectors/GD_error_dd_lambda1em4;
load ./results/error_vectors/SGD_error_dd_lambda1em4;
load ./results/error_vectors/FISTA_error_dd_lambda1em4;
load ./results/error_vectors/PDHG_error_dd_lambda1em4;
load ./results/error_vectors/SPDHG_error_dd_lambda1em4;

% Choose index
GD.plotIndex = 3;
SGD.plotIndex = 3;
FISTA.plotIndex = 2;
PDHG.plotIndex = 3;
SPDHG.plotIndex = 2;
% Plot
figure();
semilogy(0:GD.nIter{GD.plotIndex}, GD_error_dd{GD.plotIndex}, '-or', 'Linewidth', 1.5);
hold on;       
semilogy(0:SGD.nIter{SGD.plotIndex}, SGD_error_dd{SGD.plotIndex}, 'Color', 'g', 'Linewidth', 1.5);
semilogy(0:FISTA.nIter{FISTA.plotIndex}, FISTA_error_dd{FISTA.plotIndex}, 'Color', 'b', 'Linewidth', 1.5);
semilogy(0:PDHG.nIter{PDHG.plotIndex}, PDHG_error_dd{PDHG.plotIndex}, 'Color', 'm', 'Linewidth', 1.5);
semilogy(0:SPDHG.nIter{SPDHG.plotIndex}, SPDHG_error_dd{SPDHG.plotIndex}, 'Color', 'c', 'Linewidth', 1.5);
legend('FB', 'S-FB', 'AFB', 'PDHG', 'S-PDHG');
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 20 3.5 6]);
xlabel('iter/epoch');
ylabel('dual error');
set(gca,'FontSize',15);
%saveas(gcf, './figures_paper/RD09_dual_error', 'epsc');

% Relative error
figure();
semilogy(1:GD.nIter{GD.plotIndex}, (GD_error_dd{GD.plotIndex}(1:end-1)-GD_error_dd{GD.plotIndex}(2:end))./GD_error_dd{GD.plotIndex}(2:end), '-or', 'Linewidth', 1.5)
hold on;       
semilogy(1:SGD.nIter{SGD.plotIndex}, (SGD_error_dd{SGD.plotIndex}(1:end-1)-SGD_error_dd{SGD.plotIndex}(2:end))./SGD_error_dd{SGD.plotIndex}(2:end), 'Color', 'g', 'Linewidth', 1.5)
semilogy(1:FISTA.nIter{FISTA.plotIndex}, (FISTA_error_dd{FISTA.plotIndex}(1:end-1)-FISTA_error_dd{FISTA.plotIndex}(2:end))./FISTA_error_dd{FISTA.plotIndex}(2:end), 'Color', 'b', 'Linewidth', 1.5)
semilogy(1:PDHG.nIter{PDHG.plotIndex}, (PDHG_error_dd{PDHG.plotIndex}(1:end-1)-PDHG_error_dd{PDHG.plotIndex}(2:end))./PDHG_error_dd{PDHG.plotIndex}(2:end), 'Color', 'm', 'Linewidth', 1.5)
semilogy(1:SPDHG.nIter{SPDHG.plotIndex}, (SPDHG_error_dd{SPDHG.plotIndex}(1:end-1)-SPDHG_error_dd{SPDHG.plotIndex}(2:end))./SPDHG_error_dd{SPDHG.plotIndex}(2:end), 'Color', 'c', 'Linewidth', 1.5)
legend('FB', 'S-FB', 'AFB', 'PDHG', 'S-PDHG');
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 20 1e-4 1e-1]);
xlabel('iter/epoch');
ylabel('relative error');
set(gca,'FontSize',15);
%saveas(gcf, './figures_paper/RD09_relative_dual_error', 'epsc');


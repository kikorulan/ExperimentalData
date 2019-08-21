% Read data from files
%cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/RD09_finger2_doubleRes_subsampled;
cd /scratch0/NOT_BACKED_UP/frullan/ExperimentalData/RD09_finger2_doubleRes_subsampled;

clear all;
close all;

% Functions
NORM = 2.5e17;
lambda = 1e-4;
[TV, D, DTV] = TVOperators(3, 'none');
norm_distance = @(x, y) sum((x(:) - y(:)).*(x(:) - y(:)));
obj_function = @(y0, y, u0) 0.5*norm_distance(y0, NORM*y) + 0.5*lambda*TV(u0);

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
%set(gca,'FontSize', 15);
yt = get(h,'YTick');
exponent = 2;
h.Ruler.Exponent = -2;
for i = 1:length(h.TickLabels)
    h.TickLabels{i} = sprintf('%2.1f',(10^exponent)*yt(i));
end
%saveas(gcf, './figures_paper/RD09_forwardSignal', 'epsc');
%saveas(gcf, './figures_paper/RD09_forwardSignal_LF', 'epsc');

% Load Adjoint Pressure
pressure_adjoint = importdata('./input_data/pixelPressure_adjoint_3600sensors.dat', ' ', 0);
pressure_adjoint = matrix2cube(pressure_adjoint, Nz);
h = plot_projection_compact(pressure_adjoint, dx, false);
%saveas(gcf, './figures_paper/RD09_pressureAdjoint', 'epsc');

%==========================================================================================================================================================================
%===============================                                                            ===============================================================================
%===============================                     GROUND TRUTH                           ===============================================================================
%===============================                                                            ===============================================================================
%==========================================================================================================================================================================
% Forward signal - full
time_signal = importdata(['./input_data/forwardSignal_reference_14400sensors.dat'], ' ', 0);
y0_full = time_signal(2:end, :);
% Forward signal - subarray
time_signal = importdata(['./input_data/forwardSignal_reference_3600sensors.dat'], ' ', 0);
y0 = time_signal(2:end, :);
%==============================
% Full sensor array
%==============================
iter_array = [1 2 3 4 5 6 7 8 9 15 30 45 60 75 90 105 120 135 150 165 180 195 210 225 240 255 270 285 300];
subarray = zeros(120, 120);
subarray(1:2:end, 1:2:end) = 1;
GT.tau = '2e17';
GT.lambda = '1e-3';
GT_error_dd = norm_distance(y0_full, 0*y0_full);
GT_sub_error_dd = norm_distance(y0, 0*y0);
for iter = iter_array
    disp(['iter ', int2str(iter)])
    % Dual
    tSignal = importdata(['./results/ground_truth/forwardSignal_FISTA_tau', GT.tau, '_lambda', GT.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
    yi = tSignal(2:end, :);
    yi_sub = yi((subarray(:) > 0), :);
    GT_error_dd = [GT_error_dd norm_distance(y0_full, NORM*yi)];
    GT_sub_error_dd = [GT_sub_error_dd norm_distance(y0, NORM*yi_sub)];
end
% Plot dual
figure();
semilogy([0 iter_array], GT_error_dd, 'Linewidth', 1.5, 'Color', 'r')
hold on;
semilogy([0 iter_array], GT_sub_error_dd, 'Linewidth', 1.5, 'Color', 'b')
box on;
grid on;
axis([0 260 3 12]);
ax = gca;
ax.GridAlpha = 0.2;
%save ./results/error_vectors/FISTA_error_groundTruth;


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
saveas(gcf, './figures_paper/RD09_GD', 'epsc');

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
saveas(gcf, './figures_paper/RD09_S-GD', 'epsc');

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
saveas(gcf, './figures_paper/RD09_FISTA', 'epsc');

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
saveas(gcf, './figures_paper/RD09_PDHG', 'epsc');


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
saveas(gcf, './figures_paper/RD09_S-PDHG', 'png');
%==========================================================================================================================================================================
%===============================                                                            ===============================================================================
%===============================                     DISTANCE ERROR                         ===============================================================================
%===============================                                                            ===============================================================================
%==========================================================================================================================================================================
disp('******* DUAL DISTANCE ********');
%======================================================================
% Gradient descent
%======================================================================
disp('GD');
GD.tau = {'2e17', '4e17', '8e17', '1.6e18'};
GD.nIter = {30, 100, 100, 30};
%GD.nIter = {10, 100, 100, 30};
GD.lambda = '1e-4';
L = length(GD.tau);
clear GD_error_dd;
for ii = 1:L
    disp(['STEP ', int2str(ii)])
    GD_error_dd{ii} = 0.5*norm_distance(y0, 0*y0);
    for iter = 1:GD.nIter{ii}
        disp(['iter', int2str(iter)])
        % Dual error
        tSignal = importdata(['./results/forward/FB/forwardSignal_GD_tau', GD.tau{ii}, '_lambda', GD.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
        yi = tSignal(2:end, :);
        % Regularization
        pixelPressureMatrix = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau{ii}, '_lambda', GD.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
        pixelPressure = max(-1e20, matrix2cube(pixelPressureMatrix, Nz));
        % Total Error
        GD_error_dd{ii} = [GD_error_dd{ii} obj_function(y0, yi, pixelPressure)];
    end
end
% Plot dual
figure();
colors = winter(L);
for ii = 1:L
    semilogy(0:GD.nIter{ii}, GD_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 100 1 100])
legend('GD 2e17', 'GD 4e17', 'GD 8e17', 'GD 1.6e18');
save ./results/error_vectors/GD_error_dd_lambda1em4 GD_error_dd GD;
 
%======================================================================
% Stochastic Gradient descent
%======================================================================
%%  disp('S-GD');
%%  SGD.tau = {'4e17', '8e17', '1.6e18', '3.2e18'};
%%  SGD.nIter = {30, 30, 30, 30};
%%  SGD.lambda = '1e-4';
%%  SGD.batch = '100';
%%  L = length(SGD.tau);
%%  clear SGD_error_dd;
%%  for ii = 1:L
%%      disp(['STEP ', int2str(ii)])
%%      SGD_error_dd{ii} = norm_distance(y0, 0*y0);
%%      for iter = 1:SGD.nIter{ii}
%%          disp(['iter ', int2str(iter)])
%%          % Dual error
%%          tSignal = importdata(['./results/forward/S-FB/forwardSignal_S-GD_tau', SGD.tau{ii}, '_lambda', SGD.lambda, '_batch', SGD.batch, '_subepoch', int2str(iter), '.dat'], ' ', 0);
%%          yi = tSignal(2:end, :);
%%          SGD_error_dd{ii} = [SGD_error_dd{ii} norm_distance(y0, NORM*yi)];
%%      end
%%  end
%%  % Plot dual
%%  figure();
%%  colors = winter(length(SGD.tau));
%%  for ii = 1:length(SGD.tau)
%%      semilogy(0:SGD.nIter{ii}, SGD_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
%%      hold on;
%%  end
%%  box on;
%%  grid on;
%%  axis([0 30 1 10]);
%%  ax = gca;
%%  ax.GridAlpha = 0.2;
%%  legend('SGD 4e17', 'SGD 8e17', 'SGD 1.6e18', 'SGD 3.2e18');
%%  save ./results/error_vectors/SGD_error_dd_lambda1em4 SGD_error_dd SGD;

%==============================
% FISTA
%==============================
%%  disp('FISTA');
%%  FISTA.tau = {'1e17', '2e17', '4e17', '8e17', '1.6e18'};
%%  FISTA.nIter = {30, 30, 30, 30, 30};
%%  FISTA.lambda = '1e-4';
%%  L = length(FISTA.tau);
%%  clear FISTA_error_dd;
%%  for ii = 1:L
%%      disp(['STEP ', int2str(ii)])
%%      FISTA_error_dd{ii} = norm_distance(y0, 0*y0);
%%      for iter = 1:FISTA.nIter{ii}
%%          disp(['iter ', int2str(iter)])
%%          % Dual
%%          tSignal = importdata(['./results/forward/AFB/forwardSignal_FISTA_tau', FISTA.tau{ii}, '_lambda', FISTA.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
%%          yi = tSignal(2:end, :);
%%          FISTA_error_dd{ii} = [FISTA_error_dd{ii} norm_distance(y0, NORM*yi)];
%%      end
%%  end
%%  % Plot dual
%%  figure();
%%  colors = winter(length(FISTA.tau));
%%  for ii = 1:length(FISTA.tau)
%%      semilogy(0:FISTA.nIter{ii}, FISTA_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
%%      hold on;
%%  end
%%  box on;
%%  grid on;
%%  axis([0 30 1 10]);
%%  ax = gca;
%%  ax.GridAlpha = 0.2;
%%  legend('FISTA 1e17', 'FISTA 2e17', 'FISTA 4e17', 'FISTA 8e17', 'FISTA 1.6e18');
%%  save ./results/error_vectors/FISTA_error_dd_lambda1em4 FISTA_error_dd FISTA;

%======================================================================
% PDHG
%======================================================================
%%  disp('PDHG');
%%  PDHG.tau = {'2e17', '4e17', '8e17', '1.6e18', '3.2e18'};
%%  PDHG.nIter = {30, 30, 30, 30, 30};
%%  PDHG.sigma = '1';
%%  PDHG.theta = '1';
%%  PDHG.lambda = '1e-4';
%%  L = length(PDHG.tau);
%%  clear PDHG_error_dd; 
%%  for ii = 1:L
%%      disp(['STEP ', int2str(ii)])
%%      PDHG_error_dd{ii} = norm_distance(y0, 0*y0);
%%      for iter = 1:PDHG.nIter{ii}
%%          disp(['iter ', int2str(iter)])
%%          % Dual
%%          tSignal = importdata(['./results/forward/PDHG/forwardSignal_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau{ii}, '_theta', PDHG.theta, '_lambda', PDHG.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
%%          yi = tSignal(2:end, :);
%%          PDHG_error_dd{ii} = [PDHG_error_dd{ii} norm_distance(y0, NORM*yi)];
%%      end
%%  end
%%  % Plot Dual
%%  figure();
%%  colors = winter(length(PDHG.tau));
%%  for ii = 1:length(PDHG.tau)
%%      semilogy(0:PDHG.nIter{ii}, PDHG_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
%%      hold on;
%%  end
%%  box on;
%%  grid on;
%%  axis([0 30 1 10]);
%%  ax = gca;
%%  ax.GridAlpha = 0.2;
%%  legend('PDHG 2e17', 'PDHG 4e17', 'PDHG 8e17', 'PDHG 1.6e18', 'PDHG 3.2e18');
%%  save ./results/error_vectors/PDHG_error_dd_lambda1em4 PDHG_error_dd PDHG;

%==============================
% SPDHG
%==============================
%%  disp('S-PDHG');
%%  SPDHG.tau = {'2.5e16', '5e16', '1e17', '2e17', '4e17'};
%%  SPDHG.sigma = '1';
%%  SPDHG.theta = '1';
%%  SPDHG.batch = '100';
%%  SPDHG.nIter = {30, 30, 30, 30, 30};
%%  SPDHG.lambda = '1e-4';
%%  L = length(SPDHG.tau);
%%  clear SPDHG_error_dd;
%%  for ii = 1:L
%%      disp(['STEP ', int2str(ii)])
%%      SPDHG_error_dd{ii} = norm_distance(y0, 0*y0);
%%      for iter = 1:SPDHG.nIter{ii}
%%          disp(['iter ', int2str(iter)])
%%          % Dual
%%          tSignal = importdata(['./results/forward/S-PDHG/forwardSignal_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau{ii}, '_theta', SPDHG.theta, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_subepoch', int2str(iter), '.dat'], ' ', 0);
%%          yi = tSignal(2:end, :);
%%          SPDHG_error_dd{ii} = [SPDHG_error_dd{ii} norm_distance(y0, NORM*yi)];
%%      end
%%  end
%%  % Plot
%%  figure();
%%  colors = winter(length(SPDHG.tau));
%%  for ii = 1:length(SPDHG.tau)
%%      semilogy(0:SPDHG.nIter{ii}, SPDHG_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
%%      hold on;
%%  end
%%  box on;
%%  grid on;
%%  ax = gca;
%%  ax.GridAlpha = 0.2;
%%  axis([0 30 1 10]);
%%  legend('SPDHG 2.5e16', 'SPDHG 5e16', 'SPDHG 1e17', 'SPDHG 2e17', 'SPDHG 4e17');
%%  save ./results/error_vectors/SPDHG_error_dd_lambda1em4 SPDHG_error_dd SPDHG;

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


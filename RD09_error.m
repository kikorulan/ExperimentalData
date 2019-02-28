% Read data from files
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/RD09_finger2_doubleRes_subsampled;

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

%============================================================================================================================================
% FORWARD PROBLEM
%============================================================================================================================================

%==================================================
% TIME SIGNAL - REAL DATA
%==================================================
% Import data
filenameData = './input_data/forwardSignal_reference_3600sensors.dat';
timeSignal = importdata(filenameData, ' ', 0);
timeRD = timeSignal(1, :);
y0_raw = timeSignal(2:end, :);

% Import data
filenameData = './input_data/forwardSignal_reference_norm_3600sensors.dat';
timeSignal = importdata(filenameData, ' ', 0);
timeRD = timeSignal(1, :);
y0 = timeSignal(2:end, :);

%============================================================================================================================================
% ADJOINT PROBLEM
%============================================================================================================================================

%==================================================
% ADJOINT k-Wave
%==================================================
PML_size = 10;
pixelPressureMatrix = h5read('./output_data/RD09_adjoint_output_3600sensors.h5', '/p_final');
pixelPressure = max(0, pixelPressureMatrix(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size));
plot_projection(pixelPressure, dx);

%==================================================
% ADJOINT RT
%==================================================
pixelPressureMatrix = importdata('input_data/pixelPressure_adjoint_3600sensors.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

%%  %========================================================================================================================
%%  % FORWARD DATA
%%  %========================================================================================================================
%%  %==============================
%%  % Stochastic Gradient Descent
%%  %==============================
%%  %SGD.tau = {'1e16', '2e16', '5e16', '1e17', '2e17', '5e17', '1e18', '2e18', '5e18'};
%%  SGD.tau = {'5e17'};
%%  SGD.lambda = {'1e-4', '2e-4', '4e-4', '8e-4'};
%%  %SGD.lambda = {'1e-4', '2e-4', '4e-4', '8e-4'};
%%  
%%  clear yi_18;
%%  SGD.epoch  = int2str(1);
%%  for j = 1:length(SGD.lambda)
%%      for i = 1:length(SGD.tau)
%%          tSignal = importdata(['./results/forward/S-FB/forwardSignal_S-GD_tau', SGD.tau{i}, '_lambda', SGD.lambda{j}, '_batch', SGD.batch, '_epoch', SGD.epoch, '.dat'], ' ', 0);
%%          yi_18{j} = tSignal(2:end, :);
%%          figure;
%%          imagesc(yi_18{j})
%%          t = title(['S-GD - t = ', SGD.tau{i}, ', l = ', SGD.lambda{j}, ', batch = ', SGD.batch, ', epoch = ', SGD.epoch, ' - homogeneous SS']);
%%          t.Visible = 'on'; 
%%          pause(0.1);
%%      end
%%  end
%%  
%%  
%%  %========================================================================================================================
%%  % PIXEL PRESSURE
%%  %========================================================================================================================
%%  % Load Initial Pressure
%%  u0Matrix = importdata('./input_data/reconstruction_5iterations.dat', ' ', 0);
%%  u0 = matrix2cube(u0Matrix, Nz);
%%  h = plot_projection(u0, dx);
%%  a = axes;
%%  a.Visible = 'off'; 
%%  t = title('Recon lowest error');
%%  t.Visible = 'on'; 
%%  
%%  
%%  %==============================
%%  % Gradient Descent
%%  %==============================
%%  GD = [];
%%  GD.tau    = '1e18';
%%  GD.lambda = '1e-2';
%%  GD.iter = '1';
%%  pixelPressureMatrix = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', GD.iter, '.dat'], ' ', 0);
%%  pixelPressure_GD = max(0, matrix2cube(pixelPressureMatrix, Nz));
%%  plot_projection(pixelPressure_GD, dx);
%%  a = axes;
%%  a.Visible = 'off'; 
%%  t = title(['GD - t = ', GD.tau, ', l = ', GD.lambda, ', iter = ', GD.iter, ' - homogeneous SS']);
%%  t.Visible = 'on'; 
%%  
%%  %==============================
%%  % Stochastic Gradient Descent
%%  %==============================
%%  SGD = [];
%%  SGD.epoch  = '1';
%%  SGD.batch  = '800';
%%  SGD.tau    = '4e18';
%%  SGD.lambda = '1e-2';
%%  pixelPressureMatrix = importdata(['./results/adjoint/S-FB/', SGD.batch, '/pixelPressure_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_subepoch', SGD.epoch, '.dat'], ' ', 0);
%%  pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
%%  plot_projection(pixelPressure, dx);
%%  a = axes;
%%  t = title(['S-GD - t = ', SGD.tau, ', l = ', SGD.lambda, ', batch = ', SGD.batch, ', epoch = ', SGD.epoch, ' - homogeneous SS']);
%%  a.Visible = 'off'; 
%%  t.Visible = 'on'; 
%%  %caxis([0, 0.75])
%%  
%%  %========================================================================================================================
%%  % ITERATIVE RECONSTRUCTION
%%  %========================================================================================================================
%%  
%%  iter = 5;
%%  %============================================================
%%  % PARAMETERS
%%  %============================================================
%%  % GD *************************************
%%  GD = [];
%%  GD.tau    = '5e17';
%%  GD.lambda = '2e-3';
%%  GD.iter   = int2str(iter);
%%  % S-GD ***********************************
%%  SGD = [];
%%  SGD.tau    = '5e17';
%%  SGD.lambda = '1e-5';
%%  SGD.batch  = '90';
%%  SGD.epoch  = int2str(iter);
%%  % FISTA **********************************
%%  FISTA = [];
%%  FISTA.tau    = '2e17';
%%  FISTA.lambda = '2e-3';
%%  FISTA.iter   = int2str(iter);
%%  % PDHG ***********************************
%%  PDHG = [];
%%  PDHG.sigma  = '2e-1';
%%  PDHG.tau    = '1e18';
%%  PDHG.theta  = '1';
%%  PDHG.lambda = '2e-3';
%%  PDHG.iter   = int2str(iter);
%%  % S-PDHG *********************************
%%  SPDHG = [];
%%  SPDHG.sigma  = '1';
%%  SPDHG.tau    = '5e17';
%%  SPDHG.theta  = '1';
%%  SPDHG.lambda = '2e-5';
%%  SPDHG.batch  = '90';
%%  SPDHG.epoch  = int2str(iter);
%%  
%%  
%%  %==============================
%%  % Gradient Descent
%%  %==============================
%%  %GD.tau = {'1e16', '2e16', '5e16', '1e17', '2e17', '5e17', '1e18', '2e18'};
%%  GD.tau = {'2e18'};
%%  GD.lambda = {'2e-3', '5e-3', '1e-2'};
%%  GD.iter = '5';
%%  for j = 1:length(GD.lambda)
%%      for i = 1:length(GD.tau)
%%          pixelPressureMatrix = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau{i}, '_lambda', GD.lambda{j}, '_iter', GD.iter, '.dat'], ' ', 0);
%%          pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
%%          plot_projection(pixelPressure, dx);
%%          a = axes;
%%          t = title(['GD - tau = ', GD.tau{i}, ', lambda = ', GD.lambda{j}, ', iter = ', GD.iter]);
%%          a.Visible = 'off'; 
%%          t.Visible = 'on'; 
%%      end
%%      %saveas(gcf, ['./figures/RD06_GD_tau', GD.tau{i}, '_lambda', GD.lambda{j}, '_iter', GD.iter, '.fig']);
%%  end
%%  
%%  %==============================
%%  % Stochastic Gradient Descent
%%  %==============================
%%  %SGD.tau = {'1e16', '2e16', '5e16', '1e17', '2e17', '5e17', '1e18', '2e18', '5e18'};
%%  SGD.tau = {'2e18'};
%%  SGD.lambda = {'1e-2', '2e-2'};
%%  %SGD.lambda = {'5e-5', '1e-4', '2e-4', '5e-4', '1e-3'};
%%  %SGD.lambda = {'1e-4', '2e-4', '4e-4', '8e-4'};
%%  SGD.batch = '3200';
%%  
%%  SGD.epoch  = int2str(5);
%%  for j = 1:length(SGD.lambda)
%%      for i = 1:length(SGD.tau)
%%          pixelPressureMatrix = importdata(['./results/adjoint/S-FB/', SGD.batch, '/pixelPressure_S-GD_tau', SGD.tau{i}, '_lambda', SGD.lambda{j}, '_batch', SGD.batch, '_subepoch', SGD.epoch, '.dat'], ' ', 0);
%%          pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
%%          plot_projection(pixelPressure, dx);
%%          a = axes;
%%          %colorbar();
%%          t = title(['S-GD - t = ', SGD.tau{i}, ', l = ', SGD.lambda{j}, ', batch = ', SGD.batch, ', epoch = ', SGD.epoch, ' - homogeneous SS']);
%%          a.Visible = 'off'; 
%%          t.Visible = 'on'; 
%%          pause(0.1);
%%      end
%%  end
%%  %saveas(gcf, ['./figures/Example74_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_epoch', SGD.epoch, '.fig']);
%%  
%%  %========================================================================================================================
%%  % DUAL DISTANCE
%%  %========================================================================================================================
%%  
%%  
%%  
%%  nIter = 5;
%%  %====================
%%  % Gradient descent - lambda variation
%%  %====================
%%  disp('GD');
%%  GD.tau    = {'1e18'};
%%  GD.lambda = {'2e-3', '5e-3', '1e-2'};
%%  clear GD_error_dd;
%%  for ii = 1:length(GD.lambda)
%%      disp(ii)
%%      GD_error_dd{ii} = norm_distance(y0, 0*y0);
%%      for iter = 1:nIter
%%          tSignal = importdata(['./results/forward/FB/forwardSignal_GD_tau', GD.tau{1}, '_lambda', GD.lambda{ii}, '_iter', int2str(iter), '.dat'], ' ', 0);
%%          yi = tSignal(2:end, 1:size(y0, 2));
%%          GD_error_dd{ii} = [GD_error_dd{ii} norm_distance(y0, yi)];
%%          GD_error_dd{ii}
%%      end
%%  end
%%  
%%  % Plot
%%  x_axis = 0:nIter;
%%  figure();
%%  axis([0 5 3e-17 4.5e-17])
%%  colors = winter(length(GD.lambda));
%%  for ii = 1:length(GD.lambda)
%%      hold on;
%%      semilogy(x_axis, GD_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5);
%%  end
%%  legend('GD lambda 2e-3', 'GD lambda 5e-3', 'GD lambda 1e-2');
%%  title(['Dual Distance Error', GD.tau{1}]);
%%  grid on;
%%  box on;
%%  ax = gca;
%%  ax.GridAlpha = 0.2;
%%  %saveas(gcf, './figures/RD04_GD_dd.fig');
%%  
%%  %====================
%%  % Stochastig Gradient descent
%%  %====================
%%  disp('S-GD');
%%  SGD.batch = '90';
%%  SGD.tau    =  {'2e18'};
%%  %SGD.lambda = {'1e-3', '2e-3'};
%%  %SGD.lambda = {'1e-2', '2e-2', '5e-2'};
%%  %SGD.lambda = {'1e-3', '2e-3', '5e-3', '1e-2'};
%%  %SGD.lambda = {'5e-3', '1e-2'};
%%  SGD.lambda = {'1e-4', '2e-4', '5e-4', '1e-3', '2e-3', '5e-3'};
%%  %SGD.lambda = {'2e-4', '5e-4', '1e-3', '2e-3', '5e-3', '1e-2', '2e-2', '5e-2'};
%%  %SGD.lambda = {'1e-2'};
%%  nIter = 5;
%%  clear SGD_error_dd;
%%  for ii = 1:length(SGD.lambda)
%%      disp(ii)
%%      SGD_error_dd{ii} = norm_distance(y0, 0*y0);
%%      for iter = 1:nIter
%%          tSignal = importdata(['./results/forward/S-FB/', SGD.batch, '/forwardSignal_S-GD_tau', SGD.tau{1}, '_lambda', SGD.lambda{ii}, '_batch', SGD.batch, '_epoch', int2str(iter), '.dat'], ' ', 0);
%%          yi = tSignal(2:end, 1:size(y0, 2));
%%          SGD_error_dd{ii} = [SGD_error_dd{ii} norm_distance(y0, yi)];
%%          SGD_error_dd{ii}
%%      end
%%  end
%%  % Plot
%%  x_axis = 0:nIter;
%%  h = figure();
%%  colors = winter(length(SGD.lambda));
%%  colors_summer = autumn(length(SGD.lambda));
%%  semilogy(0:nIter, GD_error_dd{2}, 'Color', 'r', 'Linewidth', 1.5);
%%  hold on;
%%  for ii = 1:length(SGD.lambda)
%%      semilogy(x_axis, SGD_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5);
%%  end
%%  axis([0 5 3e-17 4.5e-17])
%%  set(gca, 'ytick', 1e-18*(30:0.5:45));
%%  %legend('GD', 'S-GD 1e-2', 'S-GD 2e-2', 'S-GD 5e-2');
%%  %legend('GD', 'S-GD 1e-3', 'S-GD 2e-3', 'S-GD 5e-3', 'S-GD 1e-2');
%%  legend('GD', 'S-GD 1e-4', 'S-GD 2e-4', 'S-GD 5e-4', 'S-GD 1e-3', 'S-GD 2e-3', 'S-GD 5e-3');
%%  %legend('GD', 'S-GD 2e-4', 'S-GD 5e-4', 'S-GD 1e-3', 'S-GD 2e-3', 'S-GD 5e-3', 'S-GD 1e-2', 'S-GD 2e-2', 'S-GD 5e-2');
%%  title(['Dual Distance Error tau = ', SGD.tau, ', batch = ', SGD.batch]);
%%  grid on;
%%  box on;
%%  ax = gca;
%%  ax.GridAlpha = 0.2;
%%  
%%  %************************************************************************************************************************
%%  %========================================================================================================================
%%  % PRIMAL ERROR
%%  %========================================================================================================================
%%  %************************************************************************************************************************
%%  % Load Initial Pressure
%%  u0Matrix = importdata('./input_data/reconstruction_5iterations.dat', ' ', 0);
%%  u0 = matrix2cube(u0Matrix, Nz);
%%  h = plot_projection(u0, dx);
%%  
%%  
%%  %====================
%%  % Gradient descent - lambda variation
%%  %====================
%%  nIter = 5;
%%  disp('GD');
%%  GD.tau    = {'2e18'};
%%  GD.lambda = {'2e-3', '5e-3', '1e-2'};
%%  clear GD_error_pd;
%%  for ii = 1:length(GD.lambda)
%%      disp(ii)
%%      GD_error_pd{ii} = norm_distance(u0, 0*u0);
%%      for iter = 1:nIter
%%          pixelPressure = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau{1}, '_lambda', GD.lambda{ii}, '_iter', int2str(iter), '.dat'], ' ', 0);
%%          pixelPressure = max(0, matrix2cube(pixelPressure, Nz));
%%          GD_error_pd{ii} = [GD_error_pd{ii} norm_distance(u0, pixelPressure)];
%%          GD_error_pd{ii}
%%      end
%%  end
%%  
%%  % Plot
%%  x_axis = 0:nIter;
%%  figure();
%%  %axis([0 10 3.6e-17 5e-17])
%%  colors = winter(length(GD.lambda));
%%  for ii = 1:length(GD.lambda)
%%      semilogy(x_axis, GD_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5);
%%      hold on;
%%  end
%%  axis([0 nIter 1 1e3]);
%%  legend('GD lambda 2e-3', 'GD lambda 5e-3', 'GD lambda 1e-2');
%%  title(['Primal Distance Error', GD.tau{1}]);
%%  grid on;
%%  box on;
%%  ax = gca;
%%  ax.GridAlpha = 0.2;
%%  %saveas(gcf, './figures/RD04_GD_dd.fig');
%%  
%%  %====================
%%  % Stochastic Gradient descent - lambda variation
%%  %====================
%%  nIter = 5
%%  disp('S-GD');
%%  SGD.batch = '800';
%%  SGD.tau    =  {'4e18'};
%%  %SGD.lambda = {'2e-3', '5e-3', '1e-2'};
%%  %SGD.lambda = {'1e-4', '2e-4', '5e-4'};
%%  SGD.lambda = {'1e-3', '2e-3', '5e-3', '1e-2'};
%%  clear SGD_error_pd;
%%  for ii = 1:length(SGD.lambda)
%%      disp(ii)
%%      SGD_error_pd{ii} = norm_distance(u0, 0*u0);
%%      for iter = 1:nIter
%%          pixelPressure = importdata(['./results/adjoint/S-FB/', SGD.batch, '/pixelPressure_S-GD_tau', SGD.tau{1}, '_lambda', SGD.lambda{ii}, '_batch', SGD.batch, '_subepoch', int2str(iter), '.dat'], ' ', 0);
%%          pixelPressure = max(0, matrix2cube(pixelPressure, Nz));
%%          SGD_error_pd{ii} = [SGD_error_pd{ii} norm_distance(u0, pixelPressure)];
%%          SGD_error_pd{ii}
%%      end
%%  end
%%  % Plot
%%  x_axis = 0:nIter;
%%  figure();
%%  colors = winter(length(SGD.lambda));
%%  semilogy(x_axis, GD_error_pd{2}, 'Color', 'r', 'Linewidth', 1.5);
%%  hold on;
%%  for ii = 1:length(SGD.lambda)
%%      semilogy(x_axis, SGD_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5);
%%  end
%%  axis([0 nIter 1 200]);
%%  legend('GD', 'S-GD 1e-4', 'S-GD 2e-4', 'S-GD 5e-4');
%%  title(['Dual Distance Error', SGD.tau, SGD.batch]);
%%  grid on;
%%  box on;
%%  ax = gca;
%%  ax.GridAlpha = 0.2;



%========================================================================================================================
% LONG RECONSTRUCTION
%========================================================================================================================
%====================
% Gradient descent - tau variation
%====================
%%  nIter = 20;
%%  disp('GD');
%%  GD.tau    = '5e17';
%%  GD.lambda = '5e-3';
%%  clear GD_error_dd;
%%  GD_error_dd = norm_distance(y0, 0*y0);
%%  for iter = 1:nIter
%%      tSignal = importdata(['./results/forward/FB/forwardSignal_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
%%      yi = tSignal(2:end, 1:end);
%%      GD_error_dd = [GD_error_dd norm_distance(y0, yi)];
%%      GD_error_dd(end)
%%  end
%%  
%%  % Plot
%%  x_axis = 0:nIter;
%%  figure();
%%  hold on;
%%  semilogy(x_axis, GD_error_dd, 'Linewidth', 1.5);
%%  %legend('GD tau 1e16', 'GD tau 2e16', 'GD tau 5e16', 'GD tau 1e17', 'GD tau 2e17', 'GD tau 5e17', 'GD tau 1e18', 'GD tau 2e18');
%%  title(['Dual Distance Error']);
%%  grid on;
%%  box on;
%%  axis([0 20 3.5e-17 5e-17])
%%  ax = gca;
%%  ax.GridAlpha = 0.2;
%%  
%%  
%%  %==============================
%%  % Gradient Descent
%%  %==============================
%%  GD.iter = '15';
%%  pixelPressureMatrix = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', GD.iter, '.dat'], ' ', 0);
%%  pixelPressure_recon = max(0, matrix2cube(pixelPressureMatrix, Nz));
%%  plot_projection(pixelPressure_recon, dx);
%%  a = axes;
%%  t = title(['GD - tau = ', GD.tau, ', lambda = ', GD.lambda, ', iter = ', GD.iter]);
%%  a.Visible = 'off'; 
%%  t.Visible = 'on';
%%  recon_15iter = cube2matrix(pixelPressure_recon);
%%  dlmwrite('./input_data/reconstruction_15iterations.dat', recon_15iter, ' ');
%%  %saveas(gcf, ['./figures/RD06_GD_tau', GD.tau{i}, '_lambda', GD.lambda{j}, '_iter', GD.iter, '.fig']);
%%  
%%  
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  disp('GD');
%%  GD.tau    = {'5e17'};
%%  GD.lambda = {'5e-3'};
%%  clear GD_error_pd;
%%  for ii = 1:length(GD.lambda)
%%      disp(ii)
%%      GD_error_pd{ii} = norm_distance(u0, 0*u0);
%%      for iter = 1:15
%%          pixelPressure = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau{1}, '_lambda', GD.lambda{ii}, '_iter', int2str(iter), '.dat'], ' ', 0);
%%          pixelPressure = max(0, matrix2cube(pixelPressure, Nz));
%%          GD_error_pd{ii} = [GD_error_pd{ii} norm_distance(u0, pixelPressure)];
%%      end
%%  end
%%  
%%  % Plot
%%  x_axis = 0:15;
%%  figure();
%%  %axis([0 10 3.6e-17 5e-17])
%%  colors = winter(length(GD.lambda));
%%  for ii = 1:length(GD.lambda)
%%      semilogy(x_axis, GD_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5);
%%      hold on;
%%  end
%%  title(['Primal Distance Error', GD.tau{1}]);
%%  grid on;
%%  box on;
%%  ax = gca;
%%  ax.GridAlpha = 0.2;
%%  %saveas(gcf, './figures/RD04_GD_dd.fig');

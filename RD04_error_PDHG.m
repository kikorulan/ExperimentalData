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
iter = 1;
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
SGD.lambda = '17.01';
SGD.batch  = '8736';
SGD.epoch  = int2str(iter);

%==============================
% Gradient Descent
%==============================
pixelPressureMatrix = importdata(['./results/pixelPressure_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', GD.iter, '.dat'], ' ', 0);
pixelPressure_GD = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure_GD, dx);
a = axes;
a.Visible = 'off'; 
t = title(['GD - t = ', GD.tau, ', l = ', GD.lambda, ', iter = ', GD.iter, ' - homogeneous SS']);
t.Visible = 'on'; 


%==============================
% Stochastic Gradient Descent
%==============================
pixelPressureMatrix = importdata(['./results/pixelPressure_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_subepoch', SGD.epoch, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['S-GD - t = ', SGD.tau, ', l = ', SGD.lambda, ', batch = ', SGD.batch, ', epoch = ', SGD.epoch, ' - homogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 


SGD.lambda = '17';
pixelPressureMatrix = importdata(['./results/adjoint/S-FB/pixelPressure_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_epoch', SGD.epoch, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['S-GD - t = ', SGD.tau, ', l = ', SGD.lambda, ', batch = ', SGD.batch, ', epoch = ', SGD.epoch, ' - homogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 

%========================================================================================================================
% DUAL DISTANCE ERROR
%========================================================================================================================
disp('******* DUAL DISTANCE ********');
% Gradient descent
disp('GD');
GD_error_dd = norm_distance(y0, 0*y0);
tSignal = importdata(['./results/forwardSignal_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
yi = tSignal(2:end, :);
GD_error_dd = norm_distance(y0, yi);
disp(['Error GD: ', GD_error_dd])


% Stochastic Gradient descent
%%  disp('S-GD');
%%  SGD_error_dd = norm_distance(y0, 0*y0);
%%  for iter = 1:5
%%      tSignal = importdata(['./results/forward/S-FB/forwardSignal_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_epoch', int2str(iter), '.dat'], ' ', 0);
%%      yi = tSignal(2:end, :);
%%      SGD_error_dd = [SGD_error_dd norm_distance(y0, yi)];
%%  end


% Plot
x_axis = [0 1 2 3 4 5 6 7 8 9 10];
figure();
plot(x_axis, GD_error_dd, 'Color', 'r', 'Linewidth', 1.5);
hold on;
%plot(x_axis, SGD_error_dd, 'Color', 'g', 'Linewidth', 1.5);
legend('FB', 'AFB', 'PDHG');
%title('Dual Distance Error');
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
saveas(gcf, ['./figures/RD04_dd_error.fig']);
saveas(gcf, ['./figures/RD04_dd_error'], 'epsc');

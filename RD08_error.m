% Read data from files
%cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/RD08_finger2_doubleRes;
cd /scratch0/NOT_BACKED_UP/frullan/ExperimentalData/RD08_finger2_doubleRes;

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
filenameData = './input_data/forwardSignal_reference_14400sensors.dat';
timeSignal = importdata(filenameData, ' ', 0);
timeRD = timeSignal(1, :);
y0_raw = timeSignal(2:end, :);

% Import data
filenameData = './input_data/forwardSignal_reference_norm_14400sensors.dat';
timeSignal = importdata(filenameData, ' ', 0);
timeRD = timeSignal(1, :);
y0 = timeSignal(2:end, :);

%============================================================================================================================================
% ADJOINT PROBLEM
%============================================================================================================================================

% ADJOINT RT
pixelPressureMatrix = importdata('output_data/pixelPressure_adjoint_RT.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title('Adjoint RT');
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, './figures/RD08_adjointRT.fig');


%========================================================================================================================
% ITERATIVE RECONSTRUCTION
%========================================================================================================================

iter = 5;
%==============================
% Gradient Descent
%==============================
GD.tau = '1e18';
GD.lambda = '5e-3';
GD.iter = '5';
pixelPressureMatrix = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', GD.iter, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['FB - tau = ', GD.tau, ', lambda = ', GD.lambda, ', iter = ', GD.iter]);
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, ['./figures/RD08_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', GD.iter, '.fig']);

%==============================
% Stochastic Gradient Descent
%==============================
SGD.tau = '2e18';
SGD.lambda = '1e-4';
SGD.batch = '90';
SGD.epoch = '5';
pixelPressureMatrix = importdata(['./results/adjoint/S-FB/', SGD.batch, '/pixelPressure_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_epoch', SGD.epoch, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['S-FB - t = ', SGD.tau, ', l = ', SGD.lambda, ', batch = ', SGD.batch, ', epoch = ', SGD.epoch, ' - homogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, ['./figures/RD08_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_epoch', SGD.epoch, '.fig']);

%========================================================================================================================
% DUAL DISTANCE
%========================================================================================================================
nIter = 5;
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
    for iter = 1:nIter
        tSignal = importdata(['./results/forward/FB/forwardSignal_GD_tau', GD.tau{1}, '_lambda', GD.lambda{ii}, '_iter', int2str(iter), '.dat'], ' ', 0);
        yi = tSignal(2:end, 1:size(y0, 2));
        GD_error_dd{ii} = [GD_error_dd{ii} norm_distance(y0, yi)];
        GD_error_dd{ii}
    end
end

% Plot
x_axis = 0:nIter;
figure();
axis([0 5 3e-17 4.5e-17])
colors = winter(length(GD.lambda));
for ii = 1:length(GD.lambda)
    hold on;
    semilogy(x_axis, GD_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5);
end
legend('FB lambda 2e-3', 'FB lambda 5e-3', 'FB lambda 1e-2');
title(['Dual Distance Error t = ' GD.tau{1}]);
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 5 3e-17 4.5e-17])
set(gca, 'ytick', 1e-18*(30:1:45));
saveas(gcf, './figures/RD08_GD_dd.fig');

%====================
% Stochastig Gradient descent
%====================
disp('S-GD');
SGD.batch = {'90', '200', '400', '1600'};
SGD.tau    =  {'2e18'};
SGD.lambda = {'1e-4', '2e-4', '5e-4', '1e-3', '2e-3', '5e-3'};
nIter = 5;

clear SGD_error_dd;
for jj = 1:length(SGD.batch)
    for ii = 1:length(SGD.lambda)
        disp(ii)
        SGD_error_dd{jj}{ii} = norm_distance(y0, 0*y0);
        for iter = 1:nIter
            tSignal = importdata(['./results/forward/S-FB/', SGD.batch{jj}, '/forwardSignal_S-GD_tau', SGD.tau{1}, '_lambda', SGD.lambda{ii}, '_batch', SGD.batch{jj}, '_subepoch', int2str(iter), '.dat'], ' ', 0);
            yi = tSignal(2:end, 1:size(y0, 2));
            SGD_error_dd{jj}{ii} = [SGD_error_dd{jj}{ii} norm_distance(y0, yi)];
            SGD_error_dd{jj}{ii}
        end
    end
    % Plot
    x_axis = 0:nIter;
    h = figure();
    colors = winter(length(SGD.lambda));
    colors_summer = autumn(length(SGD.lambda));
    for ii = 1:length(SGD.lambda)
        hold on;
        semilogy(x_axis, SGD_error_dd{jj}{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5);
    end
    axis([0 5 3e-17 4.5e-17])
    set(gca, 'ytick', 1e-18*(30:1:45));
    legend('S-FB lambda 1e-4', 'S-FB lambda 2e-4', 'S-FB lambda 5e-4', 'S-FB lambda 1e-3', 'S-FB lambda 2e-3', 'S-FB lambda 5e-3');
    t = title(['Dual Distance Error tau = ' SGD.tau{1} ', batch = ' SGD.batch{jj}]);
    grid on;
    box on;
    ax = gca;
    ax.GridAlpha = 0.2;
    pause(.1)
    saveas(gcf, ['./figures/RD08_S-GD_dd_batch' SGD.batch{jj} '.fig']);
end

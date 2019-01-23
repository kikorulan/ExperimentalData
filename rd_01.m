% CS Operator used to compress the experimantal data, here Hadamard and scrambled Hadamard patters were used
cs_operator.class = 'fast';
cs_operator.type = 'hadamard';
% Hadamard transfrom
%hH = @(x) H*x; %using explicit H matrix
hH = @(x) FHT(x); %using block fast transform: for large N

% Load full Scrambled Hadamard PAT data from file
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/Scripts
%load scrambledHadamardReal_hairknot_again.mat
load /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/DMD_PA_datasets/scrambledHadamardReal_hairknot_new.mat

% Size of the DMD: NxN
N = 128; 
% mkdir(' ../../../data/patcs/hairknot_128x128')

% Bring the data to the expected format
cs_sensor_data = permute(data_DCsub, [2 1]); %timestep second variable
figure(3), clf; plot(mean(cs_sensor_data,1)); hold on;

% Save in PATCS format
save /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/rd_01/HAIRKNOT_128x128 cperm icperm rows cs_sensor_data

% Scale and shift data
% Construct calibration function handle for the used operator
calibrateSensorData = str2func(['calibrate' capitalize(cs_operator.class)  capitalize(cs_operator.type) 'SensorData']);
% Calibrates the Hadamard data, in particular estimates the unusable all one pattern 
cs_sensor_data = calibrateSensorData(cs_sensor_data, rows, 4);
figure(3), plot(mean(cs_sensor_data,1), 'r'); legend('Raw data', 'Calibration 4')
%cs_sensor_data = calibrateSensorData(cs_sensor_data, rows, 3);
%figure(3), plot(mean(cs_sensor_data,1), 'k'); legend('Raw data', 'Calibration 2', 'Calibration 3')
%cs_sensor_data = calibrateSensorData(cs_sensor_data, rows, 1);
%figure(3), plot(mean(cs_sensor_data,1), 'g'); hold off; legend('Raw data', 'Calibration 2', 'Calibration 3', 'Calibration 1')

% Save in PATCS format
save /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/rd_01/hairknot_128x128_cssd_fast_hadamard cs_sensor_data

% Reconstruction for full scrambled Hadamard data
[~,irows] = sort(rows);
sensor_data = [];
sensor_data(icperm,:) = hH(cs_sensor_data(irows,:));
rsd = reshape(sensor_data,[N,N,size(cs_sensor_data,2)]);
scrollView(rsd,3)

% Save as gold standard data in CS PAT format
save /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/rd_01/hairknot_128x128_sd sensor_data

% Save SGL fine
dx = 5 * 13.6 *cos(24*pi/180) * 1e-6; %[m]
dy = 5 * 13.6 * 1e-6; % [m]
dt = 20e-9; %[s]
% sound_speed = 1490;
% t0 = 80;
%saveSGL(reshape(sensor_data, N, N, []), dx, dy, dt, ['/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/rd_01/hairknot_128x128_sd_fast_hadamard.SGL'])


% Reconstruction of partial scrambled Hadamard data. Needs extra calibration
% with only the used data.
M = ceil(0.18*N^2);
[~,irows] = sort(rows);
ind_all1 = find(rows==1);
ind_half1 = find(rows==(N^2/2+1));
% Make sure that all 1 pattern is always there
if ind_all1 > M
  ind_rows = [1:M-1 ind_all1];
else
  ind_rows = 1:M;
end
if ind_half1 > M
  ind_rows(M-1) = ind_half1;
end

% Bring the data to the expected format
cs_sensor_data = permute(data_DCsub, [2 1]); %timestep second variable
figure(3), clf; plot(mean(cs_sensor_data,1)); hold on;
% Scale and shift data
% Construct calibration function handle for the used operator
calibrateSensorData = str2func(['calibrate' capitalize(cs_operator.class)  capitalize(cs_operator.type) 'SensorData']);
% Calibrates the Hadamard data, in particular estimates the unusable all one pattern 
cs_sensor_data__ = calibrateSensorData(cs_sensor_data(ind_rows,:), rows(ind_rows), 4, N^2);
figure(3), plot(mean(cs_sensor_data__,1), 'r'); legend('Raw data', 'Calibration 4')

sensor_data_ = [];
cs_sensor_data_ = zeros(size(cs_sensor_data));
cs_sensor_data_(rows(ind_rows),:) = cs_sensor_data__;
sensor_data_(icperm,:) = hH(cs_sensor_data_);
rsd_ = reshape(sensor_data_,[N,N,size(cs_sensor_data_,2)]);
scrollView(rsd_,3)

% Save SGL fine
dx = 5 * 13.6 *cos(24*pi/180) * 1e-6; %[m]
dy = 5 * 13.6 * 1e-6; % [m]
dt = 20e-9; %[s]
% sound_speed = 1490
% t0 = 35
%saveSGL(reshape(sensor_data_, N, N, []), dx, dy, dt, ['/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/rd_01/hairknot_0p18_128x128_sd_fast_hadamard.SGL'])



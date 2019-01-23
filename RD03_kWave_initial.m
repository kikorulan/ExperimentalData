% THIS FILE IS AN EXAMPLE ONLY. DO NOT MODIFY THIS FILE TO RUN EXPERIMENTS!
% MAKE YOUR OWN COPY IN WHICH YOU ADJUST PATHS, OPTIONS ETC.
%
%
%
% Copyright (C) 2016 Felix Lucka


clear all; close all; clc

dataPath = '/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/RD03_vesselF/';
cd '/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/RD03_vesselF/';
setenv LD_LIBRARY_PATH '/cs/research/medim/projects2/projects/frullan/lib/root/lib64:/cs/research/medim/projects2/projects/frullan/lib/glibc-2.27/mathvec';
%setenv LD_LIBRARY_PATH '/cs/research/medim/projects2/projects/frullan/lib/root/lib64';



%%


rng('default') % Reset rand, randn, randi to default seed

% choose a data set
ID = 'BloodVessel';
settingAlias = 'Up2Iso';

startComputationInSec(0) % start the computation later

%% all the visualization parameter

visuFL = true;


visuParaRes = [];
visuParaRes.colorMap = 'parula';
visuParaRes.smoothFL = false;
visuParaRes.scaling = 'autoLinPos';
visuParaRes.scaling = 'histCutOffPos';
visuParaRes.hist_cut_off = 1/1000;
visuParaRes.printFL = true;


visuParaIter = visuParaRes;
visuParaIter.printFL = false;

imagesPostFix = ['_mI' visuParaRes.colorMap '-hist+'];

%% get all the data

% get representation of that data
dataRep = getDataRepresentation(ID,settingAlias);
data = dataRep.data;
dataRep.data = [];
dataRep.data{1} = data;
settingPara = dataRep.setting;
settingPara.computation.smoothFL = true;
settingPara.computation.dataCast = 'gpuArray-single';
settingPara.computation.codeVersion = 'kWave';
%settingPara.computation.kWaveCodeVersion = 'CUDA'; % we use the CUDA code

% generate PAT setting struct
setting = constructSetting(settingPara);
setting = initializeSensorScheme(dataRep,setting);


visuParaRes.slices2print = floor(0.2678 * setting.Ny);
visuParaRes.dimSlice     = 3;
visuParaRes.fps = 100;

% preprocessing

preProcPara = dataRep.preProcPara;
preProcPara.rmMedianFL   = false;
preProcPara.filterFL     = false;
%preProcPara.filterPara.type   = 'highPass';
%preProcPara.filterPara.frequencies = 0.2 * 10^6;
%preProcPara.filterPara.transWidth  = 0.005;
preProcPara.deleteChannelFL       = false;
%preProcPara.deleteChannelPara.percent  = 0.01;
%preProcPara.deleteChannelPara.baseline = 200:600;

preProc = constructPreProcessing(setting,preProcPara);

%% get the data and preProcess it 

getPara = [];
getPara.directPreProcFL = true;
[f,preProc] = getData(dataRep,setting,preProc,getPara);
preProc = myRmfield(preProc,'paraData');


visuSenPara = [];
visuSenPara.visuFL     = true;
visuSenPara.colorMap   = 'cool2hot';
visuSenPara.clim       = 1 * max(abs(f(:))) * [-1,1];
visuSenPara.fps = 10;
%flyThroughData(f(:,:),setting,preProc,visuSenPara)

%% define a sub-sampling pattern

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  % generate setting in which subsampling occurs 
%%  subSampleSetting = setting;
%%  subSampleSetting.sensor.subSamplingScheme = 'uniform';
%%  subSampleSetting.sensor.Msub = 1/4;% 4, 8
%%  subSampleSetting = constructSetting(subSampleSetting);
%%  subSampleSetting = initializeSensorScheme(dataRep,subSampleSetting);
%%  % sub sample data
%%  fSS = subSampleFullData(f,subSampleSetting);
%%  
%%  % generate preprocessing by zero filling
%%  decomPara = preProcPara; 
%%  decomPara.decompressionFL = true;
%%  decomPara.decompressionPara.type = 'fillIn';
%%  preProcFill0 = constructPreProcessing(subSampleSetting,decomPara);
%%  fSS0fill = preProcFill0.decompressionPara.decompressH(fSS);
%%  % this overwrite the deleted channels by 0's (but is not efficient)
%%  fSS0fill = deleteChannels(deleteChannels(fSS0fill,setting,preProc,0),setting,preProc,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% some more parameters
close all


% get the Libschitz constants of the two settings
L = getSettingProperties('LipschitzConstant',setting);
%%  L_SS = getSettingProperties('LipschitzConstant',subSampleSetting);


invParaTmp = [];
invParaTmp.spatialBC      = '0NB';
invParaTmp.constraint          = 'positivity';
invParaTmp.algo.type = 'ProximalGradient';
invParaTmp.algo.outputFL = true;
invParaTmp.algo.monitorFL = true;
invParaTmp.algo.visuFL = false;
invParaTmp.algo.visuPara   = visuParaIter;
invParaTmp.algo.visuCommand = @(x,para)  plotReconstructedSolution(x,setting,para);
invParaTmp.algo.stopCriterion    = 'maxIter';
invParaTmp.algo.breakCriterion   = 'energyIncrease';
invParaTmp.algo.breakTimes       = 10;
invParaTmp.algo.maxIter          = 50;
invParaTmp.algo.fastGradientFL   = true;
invParaTmp.algo.restartFL        = true;
invParaTmp.algo.tauUpdateRule    = 'FISTA';
invParaTmp.algo.nue              = 1.8/L;
invParaTmp.algo.monitorFL = true;
invParaTmp.BregmanIterationsFL = false;
invParaTmp.algo.denoisingAlgo.maxIter = 250;
invParaTmp.algo.denoisingAlgo.stopTolerance = 10^-10;



postProc           = [];
postProc.nonNegCon = true;
postProc.TVdenoise = false;
postProc.threshold  = 0;
postProc.blankLayer = 0;
postProc.ccComp     = 0;
postProcTV = postProc;
postProcTV.TVdenoisePara = [];
postProcTV.TVdenoise = true;
postProcTV.TVdenoisePara.maxIter = 100;
postProcTV.TVdenoisePara.BC = '0NB';
postProcTV.TVdenoisePara.stopTolerance = 10^-6;
postProcTV.TVdenoisePara.lambda = 0.15;
postProcTV.TVdenoisePara.visuFL = false;
postProcTV.TVdenoisePara.dataCast = setting.computation.dataCast;


%% time reversal
close all

% construct inverse method
invPara = invParaTmp;
invPara.type = 'TimeReversal';
invPara.algo.fillIn = 'zeros';
invMethod = constructInvMethod(dataRep,setting,preProc,invPara);
% invert full data
[p0TR,infoAux] = getInverseSolution(f,dataRep,setting,preProc,invMethod);
disp(['computation time: ' infoAux.computationTimeStr]);

% post process with TV denoising and plot and print
[p0TRpp,ppInfo] = postProcessStaticReconstruction(p0TR,setting,postProc);
[auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,setting,preProc,invMethod);
visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
maxIntensityProjection(p0TRpp,setting,visuParaRes);
%flyThroughReconstruction(p0TRpp,setting,visuParaRes);
% post process with TV denoising and plot and print
[p0TRpp,ppInfo] = postProcessStaticReconstruction(p0TR,setting,postProcTV);
[auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,setting,preProc,invMethod);
visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
maxIntensityProjection(p0TRpp,setting,visuParaRes);
%flyThroughReconstruction(p0TRpp,setting,visuParaRes);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  % invert sub sampled data
%%  [p0TRSS,infoAux] = getInverseSolution(fSS,dataRep,subSampleSetting,preProc,invMethod);
%%  disp(['computation time: ' infoAux.computationTimeStr]);
%%  % post process with TV denoising and plot and print
%%  [p0TRSSpp,ppInfo] = postProcessStaticReconstruction(p0TRSS,subSampleSetting,postProcTV);
%%  [auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,subSampleSetting,preProc,invMethod);
%%  visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
%%  maxIntensityProjection(p0TRSSpp,subSampleSetting,visuParaRes);
%%  flyThroughReconstruction(p0TRSSpp,subSampleSetting,visuParaRes);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  % invert zeor filled sub sampled data
%%  invPara.algo.fillIn = 'zeros';
%%  invMethod = constructInvMethod(dataRep,setting,preProc,invPara);
%%  [p0TRSS0,infoAux] = getInverseSolution(fSS0fill,dataRep,subSampleSetting,preProcFill0,invMethod);
%%  disp(['computation time: ' infoAux.computationTimeStr]);
%%  % post process with  and plot and print
%%  [p0TRSS0pp,ppInfo] = postProcessStaticReconstruction(p0TRSS0,subSampleSetting,postProc);
%%  [auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,subSampleSetting,preProcFill0,invMethod);
%%  visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
%%  maxIntensityProjection(p0TRSS0pp,subSampleSetting,visuParaRes);
%%  flyThroughReconstruction(p0TRSS0pp,subSampleSetting,visuParaRes);
%%  % post process with TV denoising and plot and print
%%  [p0TRSS0pp,ppInfo] = postProcessStaticReconstruction(p0TRSS0,subSampleSetting,postProcTV);
%%  [auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,subSampleSetting,preProcFill0,invMethod);
%%  visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
%%  maxIntensityProjection(p0TRSS0pp,subSampleSetting,visuParaRes);
%%  flyThroughReconstruction(p0TRSS0pp,subSampleSetting,visuParaRes);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% back projection via adjoint
close all

% construct inverse method
invPara = invParaTmp;
invPara.type = 'BackPropagation';
invMethod = constructInvMethod(dataRep,setting,preProc,invPara);

% invert full data
[p0BP,infoBP] = getInverseSolution(f,dataRep,setting,preProc,invMethod);
disp(['computation time: ' infoBP.computationTimeStr]);
% post process with TV denoising and plot and print
[p0BPpp,ppInfo] = postProcessStaticReconstruction(p0BP,setting,postProc);
[auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,setting,preProc,invMethod);
visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
maxIntensityProjection(p0BPpp,setting,visuParaRes);
%flyThroughReconstruction(p0BPpp,setting,visuParaRes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  % invert sub sampled data
%%  [p0BPSS,infoAux] = getInverseSolution(fSS,dataRep,subSampleSetting,preProc,invMethod);
%%  disp(['computation time: ' infoAux.computationTimeStr]);
%%  % post process with TV denoising and plot and print
%%  [p0BPSSpp,ppInfo] = postProcessStaticReconstruction(p0BPSS,subSampleSetting,postProc);
%%  [auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,subSampleSetting,preProc,invMethod);
%%  visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
%%  maxIntensityProjection(p0BPSSpp,subSampleSetting,visuParaRes);
%%  flyThroughReconstruction(p0BPSSpp,subSampleSetting,visuParaRes);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Solve positivity constrained least squares problem via proximal gradient algorithm

close all

% construct inverse method
invPara = invParaTmp;
invPara.type           = 'iterativeLeastSquares';
invPara.algo.stopCriterion    = 'maxIter';
invPara.algo.maxIter          = 5;
invMethod = constructInvMethod(dataRep,setting,preProc,invPara);

% invert full data
[p0PIPos,infoPI] = getInverseSolution(f,dataRep,setting,preProc,invMethod);
disp(['computation time: ' infoPI.computationTimeStr]);
% post process with TV denoising and plot and print
[p0PIPospp,ppInfo] = postProcessStaticReconstruction(p0PIPos,setting,postProc);
[auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,setting,preProc,invMethod);
visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
maxIntensityProjection(p0PIPospp,setting,visuParaRes);
%flyThroughReconstruction(p0PIPospp,setting,visuParaRes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  % invert sub sampled data
%%  invPara.algo.nue = 1.8/L_SS;
%%  invMethod = constructInvMethod(dataRep,subSampleSetting,preProc,invPara);
%%  [p0PIPosSS,infoPISS] = getInverseSolution(fSS,dataRep,subSampleSetting,preProc,invMethod);
%%  disp(['computation time: ' infoPISS.computationTimeStr]);
%%  % post process with TV denoising and plot and print
%%  [p0PIPosSSpp,ppInfo] = postProcessStaticReconstruction(p0PIPosSS,subSampleSetting,postProc);
%%  [auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,subSampleSetting,preProc,invMethod);
%%  visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
%%  maxIntensityProjection(p0PIPosSSpp,subSampleSetting,visuParaRes);
%%  flyThroughReconstruction(p0PIPosSSpp,subSampleSetting,visuParaRes);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Solve positivity constraint TV-regularized inverse Problem via proximal gradient algorithm
close all

lambdaTVreg = 5e-3;
lambdaSSFac = 4;

% construct inverse method
invPara = invParaTmp;
invPara.type           = 'VariationalRegularization';
invPara.spatialFunctional = 'TV';
invPara.lambda = lambdaTVreg;
invMethod = constructInvMethod(dataRep,setting,preProc,invPara);

% invert full data
[p0TV,infoTV] = getInverseSolution(f,dataRep,setting,preProc,invMethod);
disp(['computation time: ' infoTV.computationTimeStr]);
% post process and plot and print
[p0TVpp,ppInfo] = postProcessStaticReconstruction(p0TV,setting,postProc);
[auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,setting,preProc,invMethod);
visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
maxIntensityProjection(p0TVpp,setting,visuParaRes);
%flyThroughReconstruction(p0TVpp,setting,visuParaRes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  % invert sub sampled data
%%  invPara.algo.nue = invParaTmp.algo.nue * L/L_SS;
%%  invPara.lambda   = lambdaSSFac * lambdaTVreg * subSampleSetting.sensor.compressionOperator.size(1)/subSampleSetting.sensor.compressionOperator.size(2);
%%  invMethod = constructInvMethod(dataRep,subSampleSetting,preProc,invPara);
%%  [p0TVSS,infoTVSS] = getInverseSolution(fSS,dataRep,subSampleSetting,preProc,invMethod);
%%  disp(['computation time: ' infoTVSS.computationTimeStr]);
%%  % post process and plot and print
%%  [p0TVSSpp,ppInfo] = postProcessStaticReconstruction(p0TVSS,subSampleSetting,postProc);
%%  [auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,subSampleSetting,preProc,invMethod);
%%  visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
%%  maxIntensityProjection(p0TVSSpp,subSampleSetting,visuParaRes);
%%  flyThroughReconstruction(p0TVSSpp,subSampleSetting,visuParaRes);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

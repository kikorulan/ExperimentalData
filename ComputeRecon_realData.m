% THIS FILE IS AN EXAMPLE ONLY. DO NOT MODIFY THIS FILE TO RUN EXPERIMENTS!
% MAKE YOUR OWN COPY IN WHICH YOU ADJUST PATHS, OPTIONS ETC.
%
%
%
% Copyright (C) 2016 Felix Lucka


% clear all; close all; clc


%%
dataPath = '/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/16beamScan/Static106umScanSteps/';
cd '/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/Scripts';
cd '/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/16beamScan/Static106umScanSteps/';
addpath(dataPath);

startComputationInSec(0) % start the computation later

rng('default') % Reset rand, randn, randi to default seed

% choose a data set
ID = '16BeamScan';
%settingAlias = 'Up1Iso44x144x132';
settingAlias = 'Up2Iso80x120x120';

LipschitzComp=false; %In case you need the Lipschitz constant for iterative recon

%% all the visualization parameter


visuFL = true;

visuSenPara = [];
visuSenPara.visuFL     = true;
visuSenPara.colorMap   = 'blue2red';
visuSenPara.fps = 50;

visuParaRes = [];
visuParaRes.type = 'maxIntensity';
visuParaRes.colorMap = 'parula';
visuParaRes.smoothFL = false;
visuParaRes.scaling = 'autoLinPos';
visuParaRes.scaling = 'histCutOffPos';
visuParaRes.hist_cut_off = 1/1000;
visuParaRes.printFL = true;
visuParaRes.slices2show  = 114;
visuParaRes.slices2print = visuParaRes.slices2show;
visuParaRes.dimSlice     = 2;
visuParaRes.animatedGifFL  = true;
visuParaRes.fps = 5;
%visuParaRes.dockedFL = true;

visuParaIter = visuParaRes;
visuParaIter.printFL = false;

imagesPostFix = ['_' visuParaRes.colorMap(1) 'h+'];


%% get all the data

dataRepID='16BeamScan';

% get representation of that data
dataRep = getDataRepresentation(dataRepID,settingAlias);
settingPara = dataRep.setting;
settingPara.computation.smoothFL = false;
settingPara.computation.dataCast = 'gpuArray-single';
% settingPara.computation.dataCast = 'cpu';
settingPara.computation.codeVersion = 'internal';
%settingPara.computation.kWaveCodeVersion = 'CUDA'; % we use the CUDA code
%settingPara.computation.outputFL = true;
% generate PAT setting struct
setting = constructSetting(settingPara);
setting = initializeSensorScheme(dataRep,setting);




% preprocessing

preProcPara = dataRep.preProcPara;
preProcPara.correctBaselineFL   = false;
preProcPara.filterFL            = false;
preProcPara.filterPara.type     = 'bandPass';
%preProcPara.filterPara.frequencies = [0.5,20] * 10^6;z
%preProcPara.filterPara.transWidth  = 0.1;
preProcPara.deleteChannelFL       = false;
%preProcPara.deleteChannelPara.percent  = 0.01;
%preProcPara.deleteChannelPara.baseline = 3:15;

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




%% define a sub-sampling pattern

% generate setting in which subsampling occurs 
subSampleSetting = setting;
subSampleSetting.sensor.subSamplingScheme = 'uniform';
subSampleSetting.sensor.Msub = 1/4;% 4, 8
subSampleSetting = constructSetting(subSampleSetting);
subSampleSetting = initializeSensorScheme(dataRep,subSampleSetting);
% sub sample data
fSS = subSampleFullData(f,subSampleSetting);


% visualize and print the sub-sampling pattern
% if(visuFL)
%     [auxPath, auxFN] = getPathAndFN('sensor',1,dataRep,subSampleSetting,[],[]);
%     subSamVisuPara = [];
%     subSamVisuPara.printFL = true;
%     subSamVisuPara.fileName = [auxPath auxFN];
%     visualizeSamplingPattern(subSampleSetting,subSamVisuPara);
% end



% get the Libschitz constants of the two settings
% lipPara.tol = 10^-3;
% L    = getSettingProperties('LipschitzConstant',setting,lipPara);
% L_SS = getSettingProperties('LipschitzConstant',subSampleSetting,lipPara);

% return

%% Define forward and adjoint operators 

subSampleSettingData=subSampleSetting;
subSampleSettingData.computation.smoothFL = false;


fs = 60e6;
fc1 = 0.5e6;
fc2 = 15e6;
[c1,c2] = butter(3,[fc1/fs fc2/fs]);

F = @(data) filtfiltAlong(c1,c2,double(data),2);


Adata     = @(x) kWaveForwardOperator(x, subSampleSettingData);
A     = @(x) kWaveForwardOperator(x, subSampleSetting);
Aadj = @(x) kWaveAdjointOperator(x, subSampleSetting);


subSampleSetting.sensor.noise.sigma = 0.001;
%Start loop over batches

recSize=subSampleSetting.Nxyz;

% return




%% Lipschitz constants
% close all

if(LipschitzComp)

    % get the Libschitz constants of the two settings
    lipPara.tol = 10^-3;
    L    = getSettingProperties('LipschitzConstant',setting,lipPara);
    L_SS = getSettingProperties('LipschitzConstant',subSampleSetting,lipPara);

else
    L=1; %placeholder
end


%% some more parameters

invParaTmp = [];
invParaTmp.spatialBC      = 'NB';
invParaTmp.constraint          = 'positivity';
invParaTmp.algo.type = 'ProximalGradient';
invParaTmp.algo.outputFL = true;
invParaTmp.algo.monitorFL = true;
invParaTmp.algo.visuFL = visuFL;
invParaTmp.algo.visuPara   = visuParaIter;
invParaTmp.algo.visuCommand = @(x,para)  plotReconstructedSolution(x,setting,para);
invParaTmp.algo.stopCriterion    = 'maxIter';
invParaTmp.algo.breakCriterion   = 'energyIncrease';
invParaTmp.algo.breakTimes       = 10;
invParaTmp.algo.maxIter          = 20;
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
postProc.blankLayer = 0;%11;
postProc.ccComp     = 0;
postProcTV = postProc;
postProcTV.TVdenoisePara = [];
postProcTV.TVdenoise = false;
postProcTV.TVdenoisePara.maxIter = 100;
postProcTV.TVdenoisePara.BC = '0NB';
postProcTV.TVdenoisePara.stopTolerance = 10^-6;
postProcTV.TVdenoisePara.lambda = 0.1;
postProcTV.TVdenoisePara.visuFL = false;
postProcTV.TVdenoisePara.dataCast = setting.computation.dataCast;


%return
%% time reversal
%close all

% construct inverse method
invPara = invParaTmp;
invPara.type = 'TimeReversal';
invPara.algo.fillIn = 'zeros';
invMethod = constructInvMethod(dataRep,setting,preProc,invPara);
% invert full data
[p0TR,infoAux] = getInverseSolution(f,dataRep,setting,preProc,invMethod);
disp(['computation time: ' infoAux.computationTimeStr]);

% post process and plot and print
[p0TRpp,ppInfo] = postProcessStaticReconstruction(p0TR,setting,postProcTV);
[auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,setting,preProc,invMethod);
visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
plotReconstructedSolution(p0TRpp,setting,visuParaRes);


% invert sub sampled data
[p0TRSS,infoAux] = getInverseSolution(fSS,dataRep,subSampleSetting,preProc,invMethod);
disp(['computation time: ' infoAux.computationTimeStr]);
% post process and plot and print
[p0TRSSpp,ppInfo] = postProcessStaticReconstruction(p0TRSS,subSampleSetting,postProcTV);
[auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,subSampleSetting,preProc,invMethod);
visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
plotReconstructedSolution(p0TRSSpp,subSampleSetting,visuParaRes);


cd '/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/Scripts';
return
%% back projection via adjoint
%close all

% construct inverse method
invPara = invParaTmp;
invPara.type = 'BackPropagation';
invMethod = constructInvMethod(dataRep,setting,preProc,invPara);

% % invert full data
[p0BP,infoBP] = getInverseSolution(f,dataRep,setting,preProc,invMethod);
disp(['computation time: ' infoBP.computationTimeStr]);
% post process  and plot and print
[p0BPpp,ppInfo] = postProcessStaticReconstruction(p0BP,setting,postProcTV);
[auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,setting,preProc,invMethod);
visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
plotReconstructedSolution(p0BPpp,setting,visuParaRes);


% invert sub sampled data
[p0BPSS,infoAux] = getInverseSolution(fSS,dataRep,subSampleSetting,preProc,invMethod);
disp(['computation time: ' infoAux.computationTimeStr]);
% post process with TV denoising and plot and print
[p0BPSSpp,ppInfo] = postProcessStaticReconstruction(p0BPSS,subSampleSetting,postProcTV);
[auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,subSampleSetting,preProc,invMethod);
visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
plotReconstructedSolution(p0BPSSpp,subSampleSetting,visuParaRes);

%% Solve positivity constrained least squares problem via proximal gradient algorithm
%close all



% construct inverse method
invPara = invParaTmp;
invPara.type           = 'iterativeLeastSquares';
invPara.algo.stopCriterion    = 'maxIter';
invPara.algo.maxIter          = 5;
invMethod = constructInvMethod(dataRep,setting,preProc,invPara);

% invert full data
[p0PIPos,infoRec] = getInverseSolution(f,dataRep,setting,preProc,invMethod);
disp(['computation time: ' infoRec.computationTimeStr]);
% post process with TV denoising and plot and print
[p0PIPospp,ppInfo] = postProcessStaticReconstruction(p0PIPos,setting,postProcTV);
[auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,setting,preProc,invMethod);
visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
plotReconstructedSolution(p0PIPospp,setting,visuParaRes);


% invert sub sampled data
invPara.algo.nue = 1.8/L_SS;
invMethod = constructInvMethod(dataRep,subSampleSetting,preProc,invPara);
[p0PIPosSS,infoRecSS] = getInverseSolution(fSS,dataRep,subSampleSetting,preProc,invMethod);
disp(['computation time: ' infoRecSS.computationTimeStr]);
% post process with TV denoising and plot and print
[p0PIPosSSpp,ppInfo] = postProcessStaticReconstruction(p0PIPosSS,subSampleSetting,postProcTV);
[auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,subSampleSetting,preProc,invMethod);
visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
plotReconstructedSolution(p0PIPosSSpp,subSampleSetting,visuParaRes);


%% Solve positivity constraint TV-regularized inverse Problem via proximal gradient algorithm
%close all

lambdaTVreg = 5e-5;
lambdaSSFac = 4;



invParaTmp.algo.maxIter = 50

% construct inverse method
invPara = invParaTmp;
invPara.type           = 'VariationalRegularization';
invPara.spatialFunctional = 'TV';
invPara.lambda = lambdaTVreg;
%%% activate this to return all iterates of optimization
%invPara.algo.returnIteratesFL;
%%%
invMethod = constructInvMethod(dataRep,setting,preProc,invPara);

% invert full data
% [p0TV,infoTV] = getInverseSolution(f,dataRep,setting,preProc,invMethod);
% disp(['computation time: ' infoTV.computationTimeStr]);
% % post process and plot and print
% [p0TVpp,ppInfo] = postProcessStaticReconstruction(p0TV,setting,postProc);
% [auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,setting,preProc,invMethod);
% visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
% plotReconstructedSolution(p0TVpp,setting,visuParaRes);

% invert sub sampled data
invPara.algo.nue = invParaTmp.algo.nue * L/L_SS;
invPara.lambda   = lambdaSSFac * lambdaTVreg * subSampleSetting.sensor.compressionOperator.size(1)/subSampleSetting.sensor.compressionOperator.size(2);
invMethod = constructInvMethod(dataRep,subSampleSetting,preProc,invPara);
[p0TVSS,infoTVSS] = getInverseSolution(fSS,dataRep,subSampleSetting,preProc,invMethod);
disp(['computation time: ' infoTVSS.computationTimeStr]);
% post process and plot and print
[p0TVSSpp,ppInfo] = postProcessStaticReconstruction(p0TVSS,subSampleSetting,postProc);
[auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,subSampleSetting,preProc,invMethod);
visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
plotReconstructedSolution(p0TVSSpp,subSampleSetting,visuParaRes);

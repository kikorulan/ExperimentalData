% THIS FILE IS AN EXAMPLE ONLY. DO NOT MODIFY THIS FILE TO RUN EXPERIMENTS!
% MAKE YOUR OWN COPY IN WHICH YOU ADJUST PATHS, OPTIONS ETC.
%
%
%
% Copyright (C) 2016 Felix Lucka


clear all; close all; clc


%%

dataPath = '/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/RD04_palm/';
cd '/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/RD04_palm/';
startComputationInSec(0) % start the computation later

rng('default') % Reset rand, randn, randi to default seed

% choose a data set
ID = 'Palm1PrProc';
settingAlias = 'Up2IsoX100';


%% all the visualization parameter

visuFL = true;


visuParaRes = [];
visuParaRes.colorMap = 'gray';
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
settingPara = dataRep.setting;
settingPara.computation.smoothFL = false;
settingPara.computation.dataCast = 'gpuArray-single';
settingPara.computation.codeVersion = 'internal';
%settingPara.computation.kWaveCodeVersion = 'CUDA'; % we use the CUDA code

% generate PAT setting struct
setting = constructSetting(settingPara);
setting = initializeSensorScheme(dataRep,setting);


visuParaRes.slices2print = 60;
visuParaRes.dimSlice     = 3;
visuParaRes.fps = 100;

% preprocessing

preProcPara = dataRep.preProcPara;
preProcPara.rmMedianFL   = false;
preProcPara.filterFL     = false;
% preProcPara.filterPara.type   = 'bandPass';
% preProcPara.filterPara.frequencies = [0.5,20] * 10^6;
% preProcPara.filterPara.transWidth  = 0.1;
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
%flyThroughData(f(:,:),setting,visuSenPara)

%% define a sub-sampling pattern

% generate setting in which subsampling occurs 
subSampleSetting = setting;
%subSampleSetting.sensor.scannerType = 'multiBeam';
%subSampleSetting.sensor.nBeams = 8;
%subSampleSetting.sensor.interleaveFactor = 3;
subSampleSetting.sensor.subSamplingScheme = 'uniform';
subSampleSetting.sensor.Msub = 1/4;% 4, 8
subSampleSetting = constructSetting(subSampleSetting);
subSampleSetting = initializeSensorScheme(dataRep,subSampleSetting);
% sub sample data
fSS = subSampleFullData(f,subSampleSetting);

% generate preprocessing by zero filling
decomPara = preProcPara; 
decomPara.decompressionFL = true;
decomPara.decompressionPara.type = 'fillIn';
preProcFill0 = constructPreProcessing(subSampleSetting,decomPara);
fSS0fill = preProcFill0.decompressionPara.decompressH(fSS);
% this overwrite the deleted channels by 0's (but is not efficient)
fSS0fill = deleteChannels(deleteChannels(fSS0fill,setting,preProc,0),setting,preProc,1);


% visualize and print the sub-sampling pattern
if(visuFL)
    [auxPath, auxFN] = getPathAndFN('sensor',1,dataRep,subSampleSetting,[],[]);
    subSamVisuPara = [];
    subSamVisuPara.printFL = true;
    subSamVisuPara.fileName = [auxPath auxFN];
    visualizeSamplingPattern(subSampleSetting,subSamVisuPara);
end

%% some more parameters
close all


% get the Libschitz constants of the two settings
lipPara.tol = 10^-3;
L = getSettingProperties('LipschitzConstant',setting,lipPara);
L_SS = getSettingProperties('LipschitzConstant',subSampleSetting,lipPara);

invParaTmp = [];
invParaTmp.spatialBC      = '0NB';
invParaTmp.constraint          = 'positivity';
invParaTmp.algo.type = 'ProximalGradient';
invParaTmp.algo.outputFL = true;
invParaTmp.algo.monitorFL = true;
invParaTmp.algo.visuFL = true;
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

% post process and plot and print
[p0TRpp,ppInfo] = postProcessStaticReconstruction(p0TR,setting,postProcTV);
[auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,setting,preProc,invMethod);
visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
maxIntensityProjection(p0TRpp,setting,visuParaRes);
flyThroughReconstruction(p0TRpp,setting,visuParaRes);


% invert sub sampled data
[p0TRSS,infoAux] = getInverseSolution(fSS,dataRep,subSampleSetting,preProc,invMethod);
disp(['computation time: ' infoAux.computationTimeStr]);
% post process and plot and print
[p0TRSSpp,ppInfo] = postProcessStaticReconstruction(p0TRSS,subSampleSetting,postProcTV);
[auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,subSampleSetting,preProc,invMethod);
visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
maxIntensityProjection(p0TRSSpp,subSampleSetting,visuParaRes);
flyThroughReconstruction(p0TRSSpp,subSampleSetting,visuParaRes);


%% Solve positivity constrained iterative time reversal
close all


% construct inverse method
invPara = invParaTmp;
invPara.type           = 'iterativeTimeReversal';
invPara.extension      = 'noEx';
invPara.algo.nue            = 1;
invPara.algo.fillIn = 'zeros';
invPara.algo.stopCriterion    = 'maxIter';
%invPara.algo.maxIter          = 50;
invPara.algo.fastGradientFL   = false;
invPara.algo.breakCriterion   = 'none';
invPara.algo.monitorFL = true;
invMethod = constructInvMethod(dataRep,setting,preProc,invPara);

% invert full data
[p0iTRPos,infoPI] = getInverseSolution(f,dataRep,setting,preProc,invMethod);
disp(['computation time: ' infoPI.computationTimeStr]);
% post process with TV denoising and plot and print
[p0iTRPospp,ppInfo] = postProcessStaticReconstruction(p0iTRPos,setting,postProcTV);
[auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,setting,preProc,invMethod);
visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
maxIntensityProjection(p0iTRPospp,setting,visuParaRes);
flyThroughReconstruction(p0iTRPospp,setting,visuParaRes);


% invert sub sampled data
invMethod = constructInvMethod(dataRep,subSampleSetting,preProc,invPara);
[p0iTRPosSS,infoPISS] = getInverseSolution(fSS,dataRep,subSampleSetting,preProc,invMethod);
disp(['computation time: ' infoPISS.computationTimeStr]);
% post process with TV denoising and plot and print
[p0iTRPosSSpp,ppInfo] = postProcessStaticReconstruction(p0iTRPosSS,subSampleSetting,postProcTV);
[auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,subSampleSetting,preProc,invMethod);
visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
maxIntensityProjection(p0iTRPosSSpp,subSampleSetting,visuParaRes);
flyThroughReconstruction(p0iTRPosSSpp,subSampleSetting,visuParaRes);


%% Solve positivity constrained iterative time reversal with TV proxi
close all


lambdaTVreg = 5e-5;
lambdaSSFac = 2;


% construct inverse method
invPara = invParaTmp;
invPara.type           = 'iterativeTimeReversal';
invPara.spatialFunctional = 'TV';
invPara.lambda = lambdaTVreg;
invPara.extension      = 'noEx';
invPara.algo.fillIn = 'zeros';
invPara.algo.nue            = 1;
invPara.algo.stopCriterion    = 'maxIter';
%invPara.algo.maxIter          = 50;
invPara.algo.breakCriterion   = 'none';
invPara.algo.fastGradientFL   = false;
invPara.algo.monitorFL = true;
invMethod = constructInvMethod(dataRep,setting,preProc,invPara);

% invert full data
[p0iTRPosTV,infoPI] = getInverseSolution(f,dataRep,setting,preProc,invMethod);
disp(['computation time: ' infoPI.computationTimeStr]);
% plot and print
[auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,setting,preProc,invMethod);
visuParaRes.fileName = [auxPath auxFN imagesPostFix];
maxIntensityProjection(p0iTRPosTV,setting,visuParaRes);
flyThroughReconstruction(p0iTRPosTV,setting,visuParaRes);


invPara.lambda   = lambdaSSFac * lambdaTVreg * subSampleSetting.sensor.compressionOperator.size(1)/subSampleSetting.sensor.compressionOperator.size(2);
invMethod = constructInvMethod(dataRep,subSampleSetting,preProc,invPara);
[p0iTRPosTVSS,infoTVSS] = getInverseSolution(fSS,dataRep,subSampleSetting,preProc,invMethod);
disp(['computation time: ' infoTVSS.computationTimeStr]);
% plot and print
[auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,subSampleSetting,preProc,invMethod);
visuParaRes.fileName = [auxPath auxFN imagesPostFix];
maxIntensityProjection(p0iTRPosTVSS,subSampleSetting,visuParaRes);
flyThroughReconstruction(p0iTRPosTVSS,subSampleSetting,visuParaRes);


%% back projection via adjoint
close all

% construct inverse method
invPara = invParaTmp;
invPara.type = 'BackPropagation';
invMethod = constructInvMethod(dataRep,setting,preProc,invPara);

% invert full data
[p0BP,infoBP] = getInverseSolution(f,dataRep,setting,preProc,invMethod);
disp(['computation time: ' infoBP.computationTimeStr]);
% post process  and plot and print
[p0BPpp,ppInfo] = postProcessStaticReconstruction(p0BP,setting,postProcTV);
[auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,setting,preProc,invMethod);
visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
maxIntensityProjection(p0BPpp,setting,visuParaRes);
flyThroughReconstruction(p0BPpp,setting,visuParaRes);


% invert sub sampled data
[p0BPSS,infoAux] = getInverseSolution(fSS,dataRep,subSampleSetting,preProc,invMethod);
disp(['computation time: ' infoAux.computationTimeStr]);
% post process with TV denoising and plot and print
[p0BPSSpp,ppInfo] = postProcessStaticReconstruction(p0BPSS,subSampleSetting,postProcTV);
[auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,subSampleSetting,preProc,invMethod);
visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
maxIntensityProjection(p0BPSSpp,subSampleSetting,visuParaRes);
flyThroughReconstruction(p0BPSSpp,subSampleSetting,visuParaRes);

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
[p0PIPospp,ppInfo] = postProcessStaticReconstruction(p0PIPos,setting,postProcTV);
[auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,setting,preProc,invMethod);
visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
maxIntensityProjection(p0PIPospp,setting,visuParaRes);
flyThroughReconstruction(p0PIPospp,setting,visuParaRes);


% invert sub sampled data
invPara.algo.nue = 1.8/L_SS;
invMethod = constructInvMethod(dataRep,subSampleSetting,preProc,invPara);
[p0PIPosSS,infoPISS] = getInverseSolution(fSS,dataRep,subSampleSetting,preProc,invMethod);
disp(['computation time: ' infoPISS.computationTimeStr]);
% post process with TV denoising and plot and print
[p0PIPosSSpp,ppInfo] = postProcessStaticReconstruction(p0PIPosSS,subSampleSetting,postProcTV);
[auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,subSampleSetting,preProc,invMethod);
visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
maxIntensityProjection(p0PIPosSSpp,subSampleSetting,visuParaRes);
flyThroughReconstruction(p0PIPosSSpp,subSampleSetting,visuParaRes);



%% Solve positivity constraint TV-regularized inverse Problem via proximal gradient algorithm
close all

lambdaTVreg = 5e-5;
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
flyThroughReconstruction(p0TVpp,setting,visuParaRes);


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
maxIntensityProjection(p0TVSSpp,subSampleSetting,visuParaRes);
flyThroughReconstruction(p0TVSSpp,subSampleSetting,visuParaRes);


%% Solve positivity constrained TV-Minimization via Bregman iterations
%close all
stopScriptErr

nBregIter = 10;
lambdaBregFac = 1.25;

% construct inverse method
invPara = invParaTmp;
invPara.type                = 'VariationalRegularization';
invPara.spatialFunctional   = 'TV';
invPara.BregmanIterationsFL = true;
invPara.lambda              = lambdaBregFac*lambdaTVreg*nBregIter;
invPara.breg.stopCriterion  = 'maxIter';
invPara.breg.nBregIter      = nBregIter;
invPara.breg.outputFL    = true;
invPara.breg.visuFL      = true;
invPara.breg.visuPara    = visuParaRes;
invPara.breg.visuCommand = @(x,para)  maxIntensityProjection(postProcessStaticReconstruction(x,setting,postProc),setting,para);
invPara.breg.maxIterDec       = 0;
invPara.algo.maxIter          = 20;
invPara.algo.breakTimes       = 5;
invMethod = constructInvMethod(dataRep,setting,preProc,invPara);
[auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,setting,preProc,invMethod);
invMethod.breg.visuPara.fileName = [auxPath auxFN '_KK%K_mI' visuParaRes.colorMap 'pp'];
% invert full data
[p0TVbreg,infoTVbreg] = getInverseSolution(f,dataRep,setting,preProc,invMethod);
disp(['computation time: ' infoTVbreg.computationTimeStr]);
% post process and plot and print
[p0TVbregpp,ppInfo] = postProcessStaticReconstruction(p0TVbreg,setting,postProc);
[auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,setting,preProc,invMethod);
visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
maxIntensityProjection(p0TVbregpp,setting,visuParaRes);
flyThroughReconstruction(p0TVbregpp,setting,visuParaRes);



% invert sub sampled data
invPara.algo.nue = invParaTmp.algo.nue * L/L_SS;
invPara.lambda   = lambdaBregFac*lambdaSSFac*lambdaTVreg*nBregIter * subSampleSetting.sensor.compressionOperator.size(1)/subSampleSetting.sensor.compressionOperator.size(2);
invPara.breg.visuCommand = @(x,para)  maxIntensityProjection(postProcessStaticReconstruction(x,subSampleSetting,postProc),subSampleSetting,para);
invMethod = constructInvMethod(dataRep,subSampleSetting,preProc,invPara);
[auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,subSampleSetting,preProc,invMethod);
invMethod.breg.visuPara.fileName = [auxPath auxFN '_KK%K_mI' visuParaRes.colorMap 'pp'];
[p0TVbregSS,infoTVbregSS] = getInverseSolution(fSS,dataRep,subSampleSetting,preProc,invMethod);
disp(['computation time: ' infoTVbregSS.computationTimeStr]);
% post process and plot and print
[p0TVbregSSpp,ppInfo] = postProcessStaticReconstruction(p0TVbregSS,subSampleSetting,postProc);
[auxPath, auxFN] = getPathAndFN('invSolution',1,dataRep,subSampleSetting,preProc,invMethod);
visuParaRes.fileName = [auxPath auxFN '_' ppInfo.ID imagesPostFix];
maxIntensityProjection(p0TVbregSSpp,subSampleSetting,visuParaRes);
flyThroughReconstruction(p0TVbregSSpp,subSampleSetting,visuParaRes);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Main template script for examining experimental data and preparing the
% variables used in the rest of the toolbox
%
% Copyright (C) 2015 Felix Lucka


clear all; close all; clc

% include your adjusted setPath file here (will set global "dataPath" variable)
dataPath = '/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/RD04_palm/';
cd '/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/RD04_palm/';

rng('default') % Reset rand, randn, randi to default seed

%%

processPara  = [];    % will contain main parameters of the scan setting
preProcPara  = [];    % will contain data pre-processing parameters (e.g., filtering)
examineFL    = false;  % we want to run some data examination scripts

processPara.dataSet = 'Palm1PP'; % name of the data set in the look up file

% prepare
dataRep      = []; % will contain information about the data set
dataRep.data = {};
dataRep.type = 'experiment';
setting      = [];

%%

scannerType = 'standard';
ID = 'Palm1PrProc';
IDshort = 'Plm1PP';
realDataPathRel = 'input_data/';
setting.type = 'static';
setting.soundSpeed = 1540;

processPara.data2Examine = 1;

% get all SGL files in that dir
dataRep.data{1} = [realDataPathRel 'palm_t0equal_minus7.sgl'];

preProcPara.invertDataFL = false;
preProcPara.correctBaselineFL   = false;
preProcPara.cutSizeFL    = false;
preProcPara.cutSizePara.squareFL = false;
preProcPara.filterFL     = false;
preProcPara.filterPara.type   = 'bandPass';
preProcPara.tauStart = -7;
preProcPara.tauInc   = 1;
preProcPara.deleteChannelFL       = false;

processPara.upscalingMethod = 'freeInnerEven';

settingAlias = 'Up2IsoX100';
switch settingAlias
    case  'Up2Iso'
        preProcPara.tauEnd   = 80;
        processPara.depthResolution = 'isotropic';
        processPara.upscalingFac = 2;
        PML = [12,16,21];
    case  'Up2Iso200'
        preProcPara.tauEnd   = 200;
        processPara.depthResolution = 'isotropic';
        processPara.upscalingFac = 2;
        PML = [16,16,21];
    case  'Up2IsoX40'
        preProcPara.tauEnd   = 200;
        processPara.depthResolution = 'isotropic';
        processPara.upscalingFac = 2;
        processPara.Nx = 40;
        PML = [12,16,21];
    case  'Up2IsoX100'
        preProcPara.tauEnd   = 200; % 378
        processPara.depthResolution = 'isotropic';
        processPara.upscalingFac = 2;
        processPara.Nx = 100;
        PML = [14,16,21];
    otherwise
        error('invalid settingAlias')
end


% output
processPara
settingAlias

dataRep.ID = ID;
dataRep.IDshort = IDshort;

%% determine parameters of the PAT setting from the data


setting.sensor.scannerType = scannerType;
switch scannerType
    case {'standard','multiBeam'}
        switch dataRep.data{1}(end-3:end)
            case '.mat'
                switch setting.type
                    case 'static'
                        if(isfield(setting.sensor,'compressionOperator'))
                            switch setting.sensor.compressionOperator.type
                                case 'hadamard'
                                    switch setting.sensor.compressionOperator.class
                                        case 'fast'
                                            rawData = load([dataPath dataRep.data{1}]);
                                            f = rawData.(dataRep.variableName);
                                            setting.sensor.compressionOperator.rows = rawData.rows;
                                            setting.sensor.compressionOperator.cperm = rawData.cperm;
                                            setting.sensor.compressionOperator.icperm = rawData.icperm;
                                            if(preProcPara.flipDimFL)
                                                setting.Nt = size(f,1);
                                            else
                                                setting.Nt = size(f,2);
                                            end
                                        otherwise
                                            notImpErr
                                    end
                                otherwise
                                    notImpErr
                            end
                        else
                            notImpErr
                        end
                        setting.sensor.coverage = 'halfPlane';
                        setting.sensor.gridVoxDist = 1;
                        setting.sensor.sensorMask  = true(setting.Ny,setting.Nz);
                    case 'dynamic'
                        % the data is in one continous file
                        load([dataPath dataRep.data{1}],'parameter')
                        setting.Ny = parameter.Nx;
                        setting.Nz = parameter.Ny;
                        setting.dy = parameter.dx;
                        setting.dz = parameter.dy;
                        setting.Nt = parameter.Nt;
                        setting.dt = checkAndAssignStruct(preProcPara,'dt','>0',parameter.dt);
                        setting.dim = 3;
                        setting.sensor.coverage = 'halfPlane';
                        setting.sensor.subSamplingScheme = 'experiment';
                        setting.sensor.dataFile = dataRep.data{1};
                        setting.sensor.dataID = IDshort;
                        setting.sensor = constructSensor(setting,setting.sensor);
                        setting.ID = 'x'; setting.IDshort = 'x';
                        setting = initializeSensorScheme(dataRep,setting,[]);
                        setting = myRmfield(setting,{'ID','IDshort'});
                end
            case {'.SGL','.sgl'}
                [f, parameter] = loadSGL([dataPath dataRep.data{end}]);
                setting.Ny = parameter.Nx;
                switch scannerType
                    case 'multiBeam'
                        setting.Ny = setting.sensor.nBeams * setting.Ny;
                end
                setting.Nz = parameter.Ny;
                setting.dy = parameter.dx;
                setting.dz = parameter.dy;
                setting.Nt = parameter.Nt;
                setting.dt = checkAndAssignStruct(preProcPara,'dt','>0',parameter.dt);
                setting.sensor.coverage = 'halfPlane';
                setting.sensor.gridVoxDist = 1;
                setting.sensor.sensorMask  = true(setting.Ny,setting.Nz);
                setting.sensor.subSamplingScheme = 'fullData';
                %setting.sensor.ID = 'x'; setting.sensor.IDshort = 'x';
            otherwise
                error('Invalid file type')
        end
    case 'vShape'
        disp('construct vShape setting...')
        conVshapePara.visuFL = true;
        conVshapePara.printFL = true;
        setting = constructVshapeSetting(dataRep.data,setting,conVshapePara);
        setting.sensor.coverage = 'twoPlanes';
    otherwise
        notImpErr
end

setting
setting.sensor

%% try to load all data to check its sanity?


sanityCheckFL = checkAndAssignStruct(processPara,'sanityCheckFL','logical',false);

if(sanityCheckFL)
    disp('try to load all the data to check the files can be read...')
    
    para.dataCast = 'single';
    para.tEnd     = length(dataRep.data);
    preProc       = constructPreProcessing(setting,para);
    
    dataRep.data = dataRep.data;
    getPara      = [];
    getPara.directPreProcFL = true;
    f = getData(dataRep,setting,preProc,getPara);
    size(f)
    clear f
end


%% LOAD AND PREPROCESS DATA


%% load data

% Depending on the scanner type and aqucition, this is more or
% less complicated

selectDataRep = dataRep;
data2Examine  = checkAndAssignStruct(processPara,'data2Examine','i,>0',1);
switch scannerType
    case 'multiBeam'
        
    otherwise
        if(isfield(setting.sensor, 'subSamplingScheme') && strcmp(setting.sensor.subSamplingScheme,'experiment'))
            % sub-sampled data
            preProcPara.tStart = min(data2Examine);
            preProcPara.tEnd   = max(data2Examine);
        else
            % conventional data
            selectDataRep.data = dataRep.data(data2Examine);
            preProcPara.tEnd   = length(selectDataRep.data);
        end
end

ppSetting = setting; % copy of setting that will undergo pre-processing

% construct pre-processing struct
preProc   = constructPreProcessing(ppSetting,preProcPara);
preProc.determineSettingPropertiesFL = true;

% read in and pre-process data
getPara = [];
getPara.directPreProcFL = true;
[f,preProc] = getData(selectDataRep,ppSetting,preProc,getPara);
ppSetting   = overwriteFields(ppSetting,preProc.settingProperties);

% modifications to handle both static and dynamic data in same script
switch setting.type
    case 'static'
        size(f)
        f = {f};
        activeChannelMask{1} = true(size(f{1},1),1);
        activeChannelMask{1}(preProc.delChannel) = false;
    case 'dynamic'
        for iFrame = 1:length(f)
            activeChannelMask{iFrame} = true(size(f{iFrame},1),1);
            activeChannelMask{iFrame}(preProc.delChannel{iFrame}) = false;
        end
        size(f)
end

% special functions for V-shape scanner
switch scannerType
    case 'vShape'
        disp('*** read in single plane data')
        ppSetting1 = extractSinglePlaneSetting(setting,'Plane1');
        %ppSetting1.type = 'static';
        ppSetting2 = extractSinglePlaneSetting(setting,'Plane2');
        %ppSetting2.type = 'static';
        [f1,preProc1] = getData(selectDataRep,ppSetting1,preProc,getPara);
        [f2,preProc2] = getData(selectDataRep,ppSetting2,preProc,getPara);
        switch setting.type
            case 'static'
                size(f1)
                size(f2)
                f1 = {f1};
                f2 = {f2};
                preProc1.delChannel = {preProc1.delChannel};
                preProc2.delChannel = {preProc2.delChannel};
            case 'dynamic'
        end
    otherwise
end

%% FIRST, VARIOUS ASPECTS OF THE DATA CAN BE EXAMINED
close all

processPara.examineDigitizationFL   = examineFL;
processPara.examineSpectrumFL       = examineFL;
processPara.examine3DDataFL         = examineFL;
processPara.examine1DStatFL         = examineFL;
processPara.examine2DStatFL         = examineFL;
processPara.examineDelChannelFL     = examineFL;
processPara.examineSingleChannelFL  = examineFL;
processPara.examineMultiBeamDataFL  = examineFL & strcmp(setting.sensor.scannerType,'multiBeam');
processPara.examineOriginalDataFL   = checkAndAssignStruct(processPara,'examineOriginalDataFL','logical',false);

%% for patterend interrogation, we may rather examine the data in original sensor space
close all

if(processPara.examineOriginalDataFL)
    switch setting.sensor.interrogationScheme
        case 'patterned'
            switch setting.sensor.compressionOperator.type
                case 'hadamard'
                    % generate preprocessing by zero filling
                    dcom.decompressionFL = true;
                    dcom.decompressionPara.type =  'backProHadamard';
                    preProcPB = constructPreProcessing(setting,dcom);
                    % fill in zeros
                    for iFrame = 1:length(f)
                        f{iFrame} = preProcPB.decompressionPara.decompressH(f{iFrame});
                    end
            end
        otherwise
            notImpErr
    end
end

%% examine digitization of data
close all

examineDigitizationFL = checkAndAssignStruct(processPara,'examineDigitizationFL','logical',false);

if(examineDigitizationFL)
    disp('examine digitization.')
    figure();
    for iFrame = 1:length(f)
        hist(vec(f{iFrame}(activeChannelMask{iFrame},:)),2^16)
        drawnow();
        hold on
    end
    hold off
end

%% examine spectrum
close all

examineSpectrumFL = checkAndAssignStruct(processPara,'examineSpectrumFL','logical',false);

if(examineSpectrumFL)
    disp('examine spectrum.')
    compStatPara = [];
    compStatPara.frequencyAnalysisFL = true;
    for iFrame = 1:length(f)
        spectrumF{iFrame} = getfield(computeDataStatistics(f{iFrame}(activeChannelMask{iFrame},:),ppSetting,compStatPara),'freqAna');
    end
    
    figure();
    for iFrame = 1:length(f)
        %plotF{iFrame} = plot(spectrumF{iFrame}.freq,log10([spectrumF{iFrame}.psdxMean',...
        %    spectrumF{iFrame}.psdxMedian',spectrumF{iFrame}.psdxMin',spectrumF{iFrame}.psdxMax']));
        plotF{iFrame} = plot(spectrumF{iFrame}.freq,log10([spectrumF{iFrame}.psdxMean',...
            spectrumF{iFrame}.psdxMedian']));
        grid on
        set(plotF{iFrame}(1),'DisplayName','mean psdx f');
        set(plotF{iFrame}(2),'DisplayName','median psdx f');
        %set(plotF{iFrame}(3),'DisplayName','min psdx f');
        %set(plotF{iFrame}(4),'DisplayName','max psdx f');
        hold on
    end
    title('Periodogram Using FFT')
    xlabel('Frequency (Hz)')
    ylabel('Power/Frequency (dB/Hz)')
    drawnow();
    hold on
    
    switch scannerType
        case 'vShape'
            for iFrame = 1:length(f)
                spectrumF1{iFrame} = getfield(computeDataStatistics(f1{iFrame},ppSetting,compStatPara),'freqAna');
                spectrumF2{iFrame} = getfield(computeDataStatistics(f2{iFrame},ppSetting,compStatPara),'freqAna');
                
                plotF1{iFrame} = plot(spectrumF1{iFrame}.freq,log10([spectrumF1{iFrame}.psdxMean',spectrumF1{iFrame}.psdxMedian']));
                %plotF1{iFrame} = plot(spectrumF1{iFrame}.freq,log10([spectrumF1{iFrame}.psdxMean',spectrumF1{iFrame}.psdxMedian',spectrumF1{iFrame}.psdxMin',spectrumF1{iFrame}.psdxMax']));
                set(plotF1{iFrame}(1),'DisplayName','mean psdx f1');
                set(plotF1{iFrame}(2),'DisplayName','median psdx f1');
                %set(plotF1(3),'DisplayName','min psdx f1');
                %set(plotF1(4),'DisplayName','max psdx f1');
                
                plotF2{iFrame} = plot(spectrumF2{iFrame}.freq,log10([spectrumF2{iFrame}.psdxMean',spectrumF2{iFrame}.psdxMedian']));
                %plotF2{iFrame} = plot(spectrumF2{iFrame}.freq,log10([spectrumF2{iFrame}.psdxMean',spectrumF2{iFrame}.psdxMedian',spectrumF2{iFrame}.psdxMin',spectrumF2{iFrame}.psdxMax']));
                set(plotF2{iFrame}(1),'DisplayName','mean psdx f2');
                set(plotF2{iFrame}(2),'DisplayName','median psdx f2');
                %set(plotF2(3),'DisplayName','min psdx f2');
                %set(plotF2(4),'DisplayName','max psdx f2');
                hold on
            end
        otherwise
    end
    hold off
    drawnow();
end

%% fly through sensor data plot
close all

examine3DDataFL = checkAndAssignStruct(processPara,'examine3DDataFL','logical',false);

if(examine3DDataFL)
    disp('examine 3D data.')
    iFrame = 1;min(length(f),20);
    
    visuSenPara = [];
    visuSenPara.visuFL     = true;
    visuSenPara.colorMap   = 'blue2red';
    visuSenPara.clim = 1 * max(abs(vec(f{iFrame}(activeChannelMask{iFrame},:)))) * [-1,1];
    visuSenPara.threshold  = 0 * max(abs(vec(f{iFrame}(activeChannelMask{iFrame},:))));
    visuSenPara.fps = 25;
    
    ppSettingStat = extractStaticSetting(ppSetting,iFrame);
    preProcStatic = extractPreProc(preProc,iFrame);
    flyThroughData(f{iFrame},ppSettingStat,preProcStatic,visuSenPara)
end


%% look at single channels
close all

examineSingleChannelFL = checkAndAssignStruct(processPara,'examineSingleChannelFL','logical',false);

if(examineSingleChannelFL)
    disp('examine single channels.')
    iFrame = min(length(f),6);
    
    Nplots = 10;
    
    switch scannerType
        case 'vShape'
            pixelIndices = ceil(linspace(1,size(f1{iFrame},1),Nplots));
            figureButterfly = figure;
            axesButterfly = axes('Parent',figureButterfly,'YTick',0,'YGrid','on','GridLineStyle','-');
            box(axesButterfly,'on');hold(axesButterfly,'all');
            plot(f1{iFrame}(pixelIndices,:)','Parent',axesButterfly,'Marker','x')
            
            pixelIndices = ceil(linspace(1,size(f2{iFrame},1),Nplots));
            figureButterfly = figure;
            axesButterfly = axes('Parent',figureButterfly,'YTick',0,'YGrid','on','GridLineStyle','-');
            box(axesButterfly,'on');hold(axesButterfly,'all');
            plot(f2{iFrame}(pixelIndices,:)','Parent',axesButterfly,'Marker','x')
        otherwise
            if(isempty(preProc.delChannel))
                pixelIndices = ceil(linspace(1,size(f{iFrame},1),Nplots));
            else
                pixelIndices = setdiff(ceil(linspace(1,size(f{iFrame},1),Nplots)),preProc.delChannel{iFrame});
            end
            figureButterfly = figure;
            axesButterfly = axes('Parent',figureButterfly,'YTick',0,'YGrid','on','GridLineStyle','-');
            box(axesButterfly,'on');hold(axesButterfly,'all');
            plot(f{iFrame}(pixelIndices,:)','Parent',axesButterfly,'Marker','x')
    end
end


%% plot 1D statistics
close all

examine1DStatFL = checkAndAssignStruct(processPara,'examine1DStatFL','logical',false);

if(examine1DStatFL)
    disp('examine spatial statistics over time.')
    
    % global temporal power
    figure();
    for iFrame = 1:length(f)
        switch scannerType
            case 'vShape'
                plot(sqrt(sum(f1{iFrame}.^2,1))); hold on
                plot(sqrt(sum(f2{iFrame}.^2,1))); hold on
            otherwise
                plot(sqrt(sum(f{iFrame}(activeChannelMask{iFrame},:).^2,1)));  hold on
        end
    end
    hold off
    
    
    disp('examine channel statistics.')
    
    % global channel power
    figure();
    for iFrame = 1:length(f)
        switch scannerType
            case 'vShape'
                plot(sqrt(sum(f1{iFrame}.^2,2))); hold on
                plot(sqrt(sum(f2{iFrame}.^2,2))); hold on
            otherwise
                plot(sqrt(sum(f{iFrame}(activeChannelMask{iFrame},:).^2,2)));
        end
    end
    hold off
end




%% spatial mean and variance (in the baseline)
close all

examine2DStatFL = checkAndAssignStruct(processPara,'examine2DStatFL','logical',false);

if(examine2DStatFL)
    disp('examine temporal statistics over space.')
    
    iFrame = min(length(f),5);
    
    %window = preProcPara.deleteChannelPara.baseline;
    window = 1:481;
    climClip = 1;
    
    visuSenPara = [];
    visuSenPara.visuFL     = true;
    visuSenPara.colorMap   = 'blue2red';
    visuSenPara.scaling    = 'manualLin';
    
    switch scannerType
        case 'vShape'
            data3D1 = f1{iFrame};
            data3D1(preProc1.delChannel{iFrame},:) = 0;
            data3D1 = reshape(data3D1(ppSetting.sensor.reorder1,window),ppSetting.sensor.Mx1,ppSetting.sensor.My1,[]);
            data3D2 = f2{iFrame};
            data3D2(preProc2.delChannel{iFrame},:) = 0;
            data3D2 = reshape(data3D2(ppSetting.sensor.reorder2,window),ppSetting.sensor.Mx2,ppSetting.sensor.My2,[]);
            
            
            meanImg1 = mean(data3D1,3);
            meanImg2 = mean(data3D2,3);
            visuSenPara.clim = [-1,1] * climClip * max([max(meanImg1(:)),max(meanImg2(:))]);
            figure();image(data2RGB(meanImg1,visuSenPara)); axis image
            figure();image(data2RGB(meanImg2,visuSenPara)); axis image
            
            visuSenPara.colorMap   = 'parula';
            varImg1 = var(data3D1,0,3);
            varImg2 = var(data3D2,0,3);
            visuSenPara.clim = [0,climClip * max([max(varImg1(:)),max(varImg2(:))])];
            
            figure();image(data2RGB(varImg1,visuSenPara)); axis image
            figure();image(data2RGB(varImg2,visuSenPara)); axis image
            
            visuSenPara.clim = sqrt(visuSenPara.clim);
            figure();image(data2RGB(sqrt(varImg1),visuSenPara)); axis image
            figure();image(data2RGB(sqrt(varImg2),visuSenPara)); axis image
            
            clear data3D*
        otherwise
            ppSettingStat = extractStaticSetting(setting,iFrame);
            data3D = f{iFrame};
            data3D(statOrDynChoice(preProc.delChannel,iFrame),:) = 0;
            data3D = reshapeDataInto3D(data3D,ppSettingStat);
            data3D = data3D(:,:,window);
            meanImg = mean(data3D,3);
            visuSenPara.clim = [-1,1] * climClip * max(meanImg(:));
            figure();image(data2RGB(meanImg,visuSenPara)); axis image
            
            visuSenPara.colorMap   = 'parula';
            varImg = var(data3D,0,3);
            visuSenPara.clim = [0,climClip * max(varImg(:))];
            figure();image(data2RGB(varImg,visuSenPara)); axis image
            visuSenPara.clim = sqrt(visuSenPara.clim);
            figure();image(data2RGB(sqrt(varImg),visuSenPara)); axis image
            clear data3D
    end
end


%% visualize deleted channels
close all

examineDelChannelFL = checkAndAssignStruct(processPara,'examineDelChannelFL','logical',false);

if(examineDelChannelFL && preProc.deleteChannelFL)
    disp('examine deleted channel distribution')
    
    iFrame = min(length(f),1);
    
    switch scannerType
        case 'vShape'
            delChannel1 = zeros(size(f1{iFrame},1),1);
            delChannel1(preProc1.delChannel{iFrame}) = 1;
            delChannel1 = reshape(delChannel1(ppSetting.sensor.reorder1),ppSetting.sensor.Mx1,ppSetting.sensor.My1,[]);
            delChannel1RGB = data2RGB(delChannel1,[]);
            figure();image(delChannel1RGB); axis image
            
            delChannel2 = zeros(size(f2{iFrame},1),1);
            delChannel2(preProc2.delChannel{iFrame}) = 1;
            delChannel2 = reshape(delChannel2(ppSetting.sensor.reorder1),ppSetting.sensor.Mx1,ppSetting.sensor.My1,[]);
            delChannel2RGB = data2RGB(delChannel2,[]);
            figure();image(delChannel2RGB); axis image
        otherwise
            ppSettingStat = extractStaticSetting(setting,iFrame);
            delChannel = find(ppSettingStat.sensor.sensorMask);
            delChannel = delChannel(statOrDynChoice(preProc.delChannel,iFrame));
            delChannelMask = zeros(size(ppSettingStat.sensor.sensorMask));
            delChannelMask(delChannel) = 1;
            delChannelMaskRGB = data2RGB(delChannelMask,[]);
            figure();image(delChannelMaskRGB); axis image
    end
end




%% use quickRecon to have a look at the data?

processPara.quickReconFL   = checkAndAssignStruct(processPara,'quickReconFL','logical',false);

if(processPara.quickReconFL)
    iFrame = 1;
    
    quickReconPath = '~/Dropbox/software/pat-toolbox/'
    addpath(quickReconPath)
    try
        ppSettingStat = extractStaticSetting(setting,iFrame);
        f2D = reshapeDataInto3D(f{iFrame},ppSettingStat,0);
        [nx, ny, nt] = size(f2D);
        % write the file
        fid = fopen([dataPath 'tmp/tmp.sgl'], 'w', 'ieee-be');
        fwrite(fid, [nx, ny, nt, ppSetting.dy*1e6, ppSetting.dz*1e6, ppSetting.dt*1e9, reshape(f2D, 1, []), reshape([], 1, [])], 'single');
        fclose(fid);
        
        
        oldDir = cd([dataPath 'tmp']);
        quickRecon
    catch exception
        cd(oldDir)
        throw(exception)
    end
    cd(oldDir)
end


%% specific tests for multibeam data
close all

examineMultiBeamDataFL = checkAndAssignStruct(processPara,'examineMultiBeamDataFL','logical',false);

if(examineMultiBeamDataFL)
    disp('examine multi beam data.')
    fSingleBeams = extractSingleBeamData(f{1},ppSetting);
    
    figure();
    for iChannel = 1:length(fSingleBeams)
        meanSingleChannels(iChannel) = mean(abs(fSingleBeams{iChannel}(:)));
        stdSingleChannels(iChannel) = std(abs(fSingleBeams{iChannel}(:)));
        %hist(vec(fSingleBeams{iChannel}(:)),2^16)
        drawnow();
        hold on
    end
    hold off
end


%% NOW WE CREATE THE REST OF SETTING AND SAVE IT

% we reset some modifications by the pre processing
ppSetting.type   = setting.type;
ppSetting.sensor = setting.sensor;

% overwrite setting
setting = ppSetting;

switch scannerType
    case {'standard','multiBeam'}
        if(processPara.upscalingFac > 1)
            setting.sensor.gridVoxDist = processPara.upscalingFac;
            
            switch processPara.upscalingMethod
                case 'freeInner'
                    % the outer voxels of the sensor grid are real sensor locations
                    setting.Ny = processPara.upscalingFac * (setting.Ny-1)  + 1;
                    setting.Nz = processPara.upscalingFac * (setting.Nz-1)  + 1;
                case 'freeInnerEven'
                    % the outer voxels of the sensor grid are real sensor
                    % locations but round to next even number
                    setting.Ny = processPara.upscalingFac * (setting.Ny-1)  + 1;
                    setting.Nz = processPara.upscalingFac * (setting.Nz-1)  + 1;
                    setting.Ny = setting.Ny + mod(setting.Ny,2);
                    setting.Nz = setting.Nz + mod(setting.Nz,2);
                case 'freeMult'
                    % just multiply the number of voxels in each direction
                    setting.Ny = processPara.upscalingFac * setting.Ny;
                    setting.Nz = processPara.upscalingFac * setting.Nz;
                case {'linear'}
                    notImpErr
                    % enlarge by one less point if using linear interpolation
                    preProcPara.upScaleMethod = processPara.upscalingMethod;
                    setting.Ny = processPara.upscalingFac * (setting.Ny-1)  + 1;
                    setting.Nz = processPara.upscalingFac * (setting.Nz-1)  + 1;
                otherwise
                    notImpErr
                    preProcPara.upScaleMethod = processPara.upscalingMethod;
                    setting.Ny = processPara.upscalingFac * setting.Ny;
                    setting.Nz = processPara.upscalingFac * setting.Nz;
            end
            setting.dy = setting.dy/processPara.upscalingFac;
            setting.dz = setting.dz/processPara.upscalingFac;
        end
        
        switch processPara.depthResolution
            case 'temporalSampling'
                setting.Nx = setting.Nt;
                setting.dx = setting.soundSpeed * setting.dt;
            case 'isotropic'
                setting.dx = setting.dy;
                totalDepthReached = setting.Nt * setting.soundSpeed * setting.dt;
                NxDefault = ceil(totalDepthReached/setting.dx);
                NxDefault = NxDefault + mod(NxDefault,2);
                setting.Nx = checkAndAssignStruct(processPara,'Nx','i,>0',NxDefault);
            otherwise
                error('not a valid depthResolution value');
        end
end


%% pre-define PML size (can be changed later)

if(ischar(PML) && strcmp(PML,'opt'))
    disp('compute optimal PML.')
    
    PMLpara.mode = PML_mode;
    PMLpara.rep = 10;
    PMLpara.precision = 'gpuArray-single';
    setting.PML_size = optimalPMLsize([setting.Nx,setting.Ny,setting.Nz],[10,30],PMLpara);
else
    setting.PML_size = PML;
end

factor(setting.Nx+2*setting.PML_size(1))
factor(setting.Ny+2*setting.PML_size(2))
factor(setting.Nz+2*setting.PML_size(3))

%% check CFL

CFL = setting.soundSpeed * setting.dt / min([setting.dx,setting.dy,setting.dz])
if(CFL > 0.3)
    warning('CFL is larger than 0.3!!!')
end


%% save a copy of the script

resultsDir = makeDir([dataPath ID]);

scriptFilename = mfilename('fullpath');
if(~isempty(scriptFilename))
    copyfile([scriptFilename '.m'],[resultsDir '/CopyProcessDataScript_' settingAlias '.m']);
end


%% save all

close all

dataRep.setting = orderfields(setting);
dataRep.size = [setting.Nx,setting.Ny,setting.Nz];
dataRep.processPara = processPara;
dataRep.preProcPara = preProcPara;
dataRep = orderfields(dataRep);



sensor_data = f{1};
save ./input_data/sensor_data_8736 sensor_data;

disp('save everything.')
save([resultsDir '/dataRepresentation_' settingAlias '.mat'],'-struct','dataRep','-v7.3')

clear size
close all

% Script that reads a dynamic data set
%
% Copyright (C) 2014 Felix Lucka

clear all; close all; clc
dataPath = '/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/RD03_vesselF/';
cd '/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/RD03_vesselF/';

% original name = HB_breakdown_Predeath_590_2@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[10ns]_07s14m13h_07-09-15_avg1_2D_raw.SGL

rng('default') % Reset rand, randn, randi to default seed

% choose a data set
realDataPathRel = 'input_data/'

% will be used to distinguish the data from simulated one
type = 'experiment';
scannerType = 'standard';
invertedData = false;

% choose an ID that will be used for this data set.
ID = 'BloodVessel'
IDshort = 'BlVe'
settingAlias = 'Up2Iso'


c = 1420;

switch settingAlias
    case  'Up1Iso'
        PML_size = 16;
        tauEnd = 150;
        depthResolution = 'isotropic';
        upScaleFac = 1;
    case  'Up2Iso'
        PML_size = 17;
        tauEnd = 150;
        depthResolution = 'isotropic';
        upScaleFac = 2;
        upScaleMethod = 'freeInner';
     case  'Up3Iso'
        PML_size = 10;
        tauEnd = 150;
        depthResolution = 'isotropic';
        upScaleFac = 3;
        upScaleMethod = 'freeInner';
    otherwise
        error('invalid settingAlias')
end

quickReconFL = false;

%% read and assemble data

try
    oldDir = cd([dataPath realDataPathRel]);
    % get all SGL files in that dir
    data = listfiles('.','*.SGL');
    data = [realDataPathRel data{1}];
    fprintf(['try to load data from file ' data])
    
    [~, parameter] = loadSGL([dataPath data]);
catch exception
    cd(oldDir)
    throw(exception)
end
cd(oldDir)

dataRep.data{1} = data;
dataRep.type = type;

setting = rmfield(parameter,{'Nx','Ny','dx','dy'});
setting.type = 'static';
setting.sensor.coverage = 'halfPlane';
setting.sensor.scannerType = scannerType;
setting.sensor.gridVoxDist = 1;
setting.sensor.sensorMask = true(parameter.Nx,parameter.Ny);
f2D = readData(setting,dataRep,[],1);

if(invertedData)
    f2D = -f2D;
end

%% examine digitization 

figure();
hist(f2D(:),2^16)

%% filter the data

preProcPara = [];
preProcPara.invertDataFL = invertedData;
preProcPara.bpFilter = false;
preProcPara.deleteChannelFL = false;
%preProcPara.bpFrequencies = [0.005 0.5];
%preProcPara.bpFiltOrder   = 4;

if(preProcPara.bpFilter)
    [c1,c2] = butter(preProcPara.bpFiltOrder,preProcPara.bpFrequencies);
    f2D = filter(c1,c2,f2D,[],3);
end
f = reshape(f2D,[],size(f2D,3));

%% start quickRecon to determine certain settings

if(quickReconFL)
    try
        [nx, ny, nt] = size(f2D);
        % write the file
        fid = fopen([dataPath 'tmp/tmp.sgl'], 'w', 'ieee-be');
        fwrite(fid, [nx, ny, nt, parameter.dx*1e6, parameter.dy*1e6, parameter.dt*1e9, reshape(f2D, 1, []), reshape([], 1, [])], 'single');
        fclose(fid);


        oldDir = cd([dataPath 'tmp']);
        quickRecon
    catch exception
        cd(oldDir)
        throw(exception)
    end
    cd(oldDir)
end


%% construct parts of the setting

setting.sensor.coverage = 'halfPlane';
setting.sensor.scannerType = scannerType;
setting.Ny = parameter.Nx;
setting.Nz = parameter.Ny;
setting.dy = parameter.dx;
setting.dz = parameter.dy;
setting.soundSpeed = c;


%% fly through data plot
close all


visuSenPara = [];
visuSenPara.visuFL     = true;
visuSenPara.colorMap   = 'blue2red';
visuSenPara.clim = 1 * max(abs(f2D(:))) * [-1,1];
visuSenPara.threshold  = 0 * max(abs(f2D(:)));
visuSenPara.fps = 5;
%flyThroughData(f2D,setting, preProcPara, visuSenPara)


%% butterfly plots
close all

Nplots = 20;

gridVec = linspace(0,1,2*sqrt(Nplots)+1);
gridVec = gridVec(2:2:end);
[Xgrid,Ygrid] = meshgrid(gridVec,gridVec);
points = [Xgrid(:),Ygrid(:)];
pixelIndices = find(pixelHit(parameter.Nx,points));
figure();



figure1 = figure;
axes1 = axes('Parent',figure1,'YTick',0,'YGrid','on','GridLineStyle','-');
box(axes1,'on');hold(axes1,'all');
plot(f(pixelIndices,:)','Parent',axes1,'Marker','x')

%% plot temporal statistics
close all

% global temporal power
figure();plot(sqrt(sum(f.^2,1)));

% spatial mean and variance
meanImg = mean(f2D,3);
figure();image(data2RGB(meanImg,visuSenPara)); axis image 
varImg = var(f2D,0,3);
figure();image(data2RGB(varImg,visuSenPara)); axis image 
stdImg = std(f2D,0,3);
figure();image(data2RGB(stdImg,visuSenPara)); axis image 





%% define start and end point

preProcPara.tauStart = 10;
preProcPara.tauEnd   = tauEnd;





%% determine preprocessing settings

if(exist('cutSize','var'))
    preProcPara.cut2size      = true;
    setting.Ny            = cutSize(1);
    setting.Nz            = cutSize(2);
else
    preProcPara.cut2size      = false;
end
preProcPara.cut2square    = true;
if(preProcPara.cut2square)
    Nplane = min(setting.Ny,setting.Nz);
    setting.Ny = Nplane;
    setting.Nz = Nplane;
end

%% construct setting
close all

setting.Nt = length(preProcPara.tauStart:preProcPara.tauEnd);
setting.PML_size = PML_size;

setting.sensor.gridVoxDist = upScaleFac;
preProcPara.upScaleMethod = 'none';
if(upScaleFac > 1)
    switch upScaleMethod
        case 'freeInner'
            % the outer voxels of the sensor grid are real sensor locations
            setting.Ny = upScaleFac * (setting.Ny-1)  + 1;
            setting.Nz = upScaleFac * (setting.Nz-1)  + 1;
        case 'freeEven'
            % just multiply the number of voxels in each direction
            setting.Ny = upScaleFac * setting.Ny;
            setting.Nz = upScaleFac * setting.Nz;
        case {'linear'}
            % enlarge by one less point if using linear interpolation
            preProcPara.upScaleMethod = upScaleMethod;
            setting.Ny = upScaleFac * (setting.Ny-1)  + 1;
            setting.Nz = upScaleFac * (setting.Nz-1)  + 1;
        otherwise
            preProcPara.upScaleMethod = upScaleMethod;
            setting.Ny = upScaleFac * setting.Ny;
            setting.Nz = upScaleFac * setting.Nz;
    end
    setting.dy = setting.dy/upScaleFac;
    setting.dz = setting.dz/upScaleFac;
end

switch depthResolution
    case 'temporalSampling'
        setting.Nx = setting.Nt;
        setting.dx = c * setting.dt;
    case 'isotropic'
        setting.dx = setting.dy;
        totalDepthReached = setting.Nt * c * setting.dt;
        setting.Nx = ceil(totalDepthReached/setting.dx);
    otherwise
        error('not a valid depthResolution value');
end

size = [setting.Nx , setting.Ny , setting.Nz];
[factor(size(1)+2*PML_size),factor(size(2)+2*PML_size),factor(size(3)+2*PML_size)]
max([factor(size(1)+2*PML_size),factor(size(2)+2*PML_size),factor(size(3)+2*PML_size)])

%% specify the sub-sampling pattern used

isSubSampled = false;

%% save everything

setting = orderfields(setting)
resultsDir = makeDir([dataPath ID]);
save([resultsDir '/dataRepresentation_' settingAlias '.mat'],'ID','IDshort','size','data',...
    'depthResolution','setting','isSubSampled','preProcPara','type')





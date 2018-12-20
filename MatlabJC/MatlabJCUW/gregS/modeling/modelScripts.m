%% Bipolar subunits 1D spatial, no time

%parameters
%nBipolars = 10;
bipolarSpacing = 80;
bipolarRFWidth = 100;
%GCRFWidth = 500;
stimSize = 1000;
rectificationThreshold = 0.4;

%set up bipolars
bpCenters = linspace(1,stimSize,floor(stimSize/bipolarSpacing));
nBipolars = length(bpCenters);

for i=1:nBipolars;
    i
    B(i) = Subunit;
    B(i).inputGain = [];
    B(i).outputNL = @(x,T)rectify(x,rectificationThreshold,NullTransform);
    
    %make filter
    f = zeros(1,stimSize);
    startPoint = bpCenters(i)-bipolarRFWidth/2+1;
    endPoint = bpCenters(i)+bipolarRFWidth/2;
    
    %edge effects
    startPoint_corrected = max([startPoint,1]);
    endPoint_corrected = min([endPoint,stimSize]);
    
    G = gausswin(bipolarRFWidth)';
    f(startPoint_corrected:endPoint_corrected) = G(startPoint_corrected - startPoint + 1:bipolarRFWidth-(endPoint - endPoint_corrected));  
    B(i).filter = f;
end

%bipolar layer
bipolarLayer = Layer;
bipolarLayer.subunits = B;

%set up ganglion cell
RGC = Subunit;
RGC.inputGain = gausswin(nBipolars)'; %gaussian summation over bipolars
RGC.outputNL = @(x,T)x; %linear
%RGC.outputNL = @(x,T)sum(x); %integrate
RGC.filter = ones(1,stimSize); %no additional filtering step

%RGC layer
RGCLayer = Layer;
RGCLayer.subunits = RGC;

%stimulus 
s = zeros(1,stimSize);
s(400:600) = 1;

%set up model
M = FFModel;
M.layers = [bipolarLayer, RGCLayer];
M.nTrials = 1;
M.stimulus = s;
M.transform = NullTransform;

M.run;

%% Bipolar subunits 2D spatial, hexagonal grid, no time

% parameters
stimSize = 200; %square, (pixels) 
bipolarSpacing = 10;
bipolarRFWidth = 8; %standard deviation
RGCRFWidth = 111; %standard deviation

textureSigma = 20;
rectificationThreshold = 0;

%set up bipolars
% Generate hexagonal grid
gridSize = ceil((stimSize./bipolarSpacing) * 1.1);
if rem(gridSize,2) == 1, gridSize = gridSize+1; end

Rad3Over2 = sqrt(3) / 2;
[X Y] = meshgrid(1:1:gridSize);
n = size(X,1);
X = Rad3Over2 * X;
Y = Y + repmat([0 0.5],[n,n/2]);

%set spacing
X = X * bipolarSpacing;
Y = Y * bipolarSpacing;

Ind = (X-bipolarRFWidth <= stimSize & Y-bipolarRFWidth <= stimSize);
X = X(Ind);
Y = Y(Ind);

%X,Y are now coordinates of bipolar centers
bpCenters = [X, Y];
nBipolars = numel(X);

G = fspecial('gaussian',bipolarRFWidth*6,bipolarRFWidth);

%transform
T = TwoD2OneDTransform;
T.rows = stimSize;
T.cols = stimSize;

for i=1:nBipolars;
    i
    B(i) = Subunit;
    B(i).inputGain = [];
    %B(i).outputNL = @(x,T)sumAndRectify(x,rectificationThreshold,T);
    B(i).outputNL = @(x,T)sum(x); %linear summation
    %B(i).outputNL = @(x,T)rectify(x,rectificationThreshold,T);
    
    
    %make filter
    f = zeros(stimSize,stimSize);
    startPointX = ceil(bpCenters(i,1)-bipolarRFWidth*3)+1; %3s.d of each side
    endPointX = floor(bpCenters(i,1)+bipolarRFWidth*3);    
    startPointY = ceil(bpCenters(i,2)-bipolarRFWidth*3)+1; %3s.d of each side
    endPointY = floor(bpCenters(i,2)+bipolarRFWidth*3);
    
    %edge effects
    startPoint_correctedX = max([startPointX,1]);
    endPoint_correctedX = min([endPointX,stimSize]);
    startPoint_correctedY = max([startPointY,1]);
    endPoint_correctedY = min([endPointY,stimSize]);
    
    fIndX = startPoint_correctedX:endPoint_correctedX;
    fIndY = startPoint_correctedY:endPoint_correctedY;
    GIndX = startPoint_correctedX - startPointX + 1:bipolarRFWidth*6-(endPointX - endPoint_correctedX);
    GIndY = startPoint_correctedY - startPointY + 1:bipolarRFWidth*6-(endPointY - endPoint_correctedY); 
    GIndX = GIndX(1:length(fIndX));
    GIndY = GIndY(1:length(fIndY));
    f(fIndX,fIndY) = G(GIndX,GIndY); 
    B(i).filter = T.invert(f);
end

%bipolar layer
bipolarLayer = Layer;
bipolarLayer.subunits = B;
bipolarLayer.transform = T;

%set up ganglion cell
RGC = Subunit;
%gaussian summation over bipolars
RGCRFWidth_bipolarUnits = RGCRFWidth / bipolarSpacing;

RGCinputs2D = fspecial('gaussian', ceil(sqrt(nBipolars)), RGCRFWidth_bipolarUnits);
Tbp = TwoD2OneDTransform; %in units of bipolars
Tbp.rows = ceil(sqrt(nBipolars));
Tbp.cols = ceil(sqrt(nBipolars));
temp = Tbp.invert(RGCinputs2D);
RGC.inputGain = temp(1:nBipolars);
RGC.outputNL = @(x,T)x; %linear
%RGC.outputNL = @(x,T)sum(x); %integrate

RGC.filter = 1; %no additional filtering step
%RGC.filter = T.invert(ones(stimSize,stimSize)); %no additional filtering step

%RGC layer
RGCLayer = Layer;
RGCLayer.subunits = RGC;
RGCLayer.transform = NullTransform;

for i=1:20
    i
%stimulus 
s = generateTexture(stimSize,textureSigma);

%set up model
M = FFModel;
M.layers = [bipolarLayer, RGCLayer];
M.nTrials = 1;
M.stimulus = T.invert(s);

M.run;

a(i) = M.output{1};
end


%% using ModelExperiment

clear E;
modelInitParameters.stimSize = 200;
modelInitParameters.bipolarSpacing = 10;
modelInitParameters.nTrials = 1;
modelInitParameters.bipolarRFWidth = {4 8 12 16 20 24};
modelInitParameters.RGCRFWidth = 111;
 
%modelFilterParams{1}.bipolarRFWidth = {4 8 12 16 20 24};
modelFilterParams{1} = struct;
modelFilterParams{2}.RGCRFWidth = 111;

%modelNLParams{1}.outputNL = {@(x,T)sum(x);  @(x,T)sumAndRectify(x,0,T)};
%modelNLParams{1}.outputNL = {@(x,T)sum(x);  @(x,T)sumAndHill(x,.7,5.7,T)};
modelNLParams{1}.outputNL = @(x,T)sumAndHill(x,.7,6,T);

modelNLParams{2}.outputNL = @(x,T)sum(x);

modelStimParams.textureSigma = {0 4 8 12 16 20 30 40};

E = ModelExperiment;
E.modelNLParams = modelNLParams;
E.modelFilterParams = modelFilterParams;
E.modelInitParams = modelInitParameters;
E.modelStimParams = modelStimParams;
E.modelInitFunction = @(params)makeBipolarSpatialSamplingModel(params);
E.setFilterFunctions{1} = @(layer,params)setInputLayerGaussianSpatialFilter(layer,params);
E.setFilterFunctions{2} = @(layer,params)NullFilter(layer,params);
E.setNLFunctions{1} = @(layer,params)uniformNLSetter(layer,params);
E.setNLFunctions{2} = @(layer,params)uniformNLSetter(layer,params);
E.setStimFunction = @(model,params)setTextureStimulus(model,params);
E.runsPerCondition = 50;
E.buildconditionMatrix
E.run;

save('ModelExperimentResults_bipolarRF.mat', 'E');

[nModelTypes, nStims] = size(E.conditionMatrix);
micronsPerPixel = 1.8;
textureSigma = zeros(1,nStims);
Model_means = zeros(nModelTypes,nStims);
Model_errs = zeros(nModelTypes,nStims);
for i=1:nModelTypes
    for j=1:nStims
        textureSigma(j) = E.conditionMatrix{i,j}.textureSigma;
        Model_means(i,j) = mean(cell2mat(E.conditionMatrix{i,j}.result));
        Model_errs(i,j) = std(cell2mat(E.conditionMatrix{i,j}.result))./sqrt(E.runsPerCondition);
    end
end
%LinModel_means = LinModel_means./max(LinModel_means);
%LinModel_errs = LinModel_errs./max(LinModel_means);
%Model_means = Model_means./repmat(max(Model_means,[],2),1,nStims);
Model_means = Model_means./repmat(Model_means(:,nStims),1,nStims);
Model_means = Model_means';
Model_errs = Model_errs';
%NLinModel_errs = NLinModel_errs./max(NLinModel_means);

textureSigma = textureSigma*micronsPerPixel;
textureSigma = repmat(textureSigma,nModelTypes,1)';

figure;
%errorbar(textureSigma,LinModel_means,LinModel_errs,'bx-');
hold on;
errorbar(textureSigma,Model_means,Model_errs,'x-');

%% Single Cone Test

T = FFT1DTransform;

Filter = load('testFilter.mat','f');
stim = load('testStim.mat','stim');
stim = stim.stim;
Filter = Filter.f;
cone = Subunit;
cone.filter = T.invert(Filter);
cone.outputNL = @(x,T)x;
cone.inputGain = [];


L = Layer;
L.transform = T;
L.subunits = cone;

M = FFModel;
M.layers = L;
M.nTrials = 1;
M.stimulus = T.invert(stim);

%% Cone Array 
clear E;
E = ModelExperiment;
%model init params and function
E.modelInitParams.stimSize = 50;
E.modelInitParams.coneSpacing = 12;
E.modelInitParams.nTrials = 20;
E.modelInitParams.coneRFWidth = 10;
E.modelInitFunction = @(params)makeConeSpatialSamplingModel(params);
%filter params and function
load LinearFilter.mat
f = LinearFilter{1};
f_dec = decimate(f,10); %change sampleRate to 1kHz 
f_dec = f_dec(1:1000); %only 1 sec stim (for now)
Ttime = FFT1DTransform;
f_fft = Ttime.invert(f_dec);
E.modelFilterParams{1}.filter = f_fft;
E.setFilterFunctions{1} = @(layer,params)uniformFilterSetter(layer,params);
%NL params and function
E.modelNLParams{1}.outputNL =  @(x,T)addConeNoise(x,1000,T);%linear output with additive noise
E.setNLFunctions{1} = @(layer,params)uniformNLSetter(layer,params);
%stim params and function
%E.modelStimParams.barLocation = 40:60;
E.modelStimParams.barLocation = 20:30;
E.modelStimParams.barFrames = 501:600;
E.modelStimParams.nFrames = 1000;
E.modelStimParams.amplitude = 1E4; 
E.setStimFunction = @(model,params)setFlashedBarStimulus(model,params);
%finish setting up E
E.runsPerCondition = 1;
E.buildconditionMatrix;

disp('Running Model');
E.run;

%collect response (just 1 for now)
meanR = zeros(E.modelInitParams.nTrials, E.modelStimParams.nFrames);
for r = 1:E.modelInitParams.nTrials
    nCones = length(E.conditionMatrix{1}.result{1}{r});
    R = zeros(nCones,E.modelStimParams.nFrames);
    for i=1:nCones
        R(i,:) = Ttime.revert(E.conditionMatrix{1}.result{1}{r}{i});
    end
    meanR(r,:) = mean(R);
end
gMean = mean(meanR);
X = 1:E.modelStimParams.nFrames;
plot(X, gMean, 'k-');
hold on;
errR = std(meanR);
plot(X,gMean + errR,'k--');
plot(X,gMean - errR,'k--');
hold off;

















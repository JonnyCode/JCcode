function M = makeBipolarSpatialSamplingModel(params)
%parameters
stimSize = params.stimSize; %square, (pixels)
bipolarSpacing = params.bipolarSpacing;
bipolarRFWidth = params.bipolarRFWidth;
RGCRFWidth = params.RGCRFWidth;
nTrials = params.nTrials;


%set up bipolars
%Generate hexagonal grid
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

%transform
T = TwoD2OneDTransform;
T.rows = stimSize;
T.cols = stimSize;

disp(['Making ' num2str(nBipolars) ' bipolars...']);
for i=1:nBipolars;
    B(i) = Subunit;
    B(i).inputGain = [];
end
disp('done');

%bipolar layer
bipolarLayer = Layer;
bipolarLayer.subunits = B;
bipolarLayer.transform = T;

%set up ganglion cell
RGC = Subunit;

%RGC layer
RGCLayer = Layer;
RGCLayer.subunits = RGC;
RGCLayer.transform = NullTransform;

%RGC input gain gaussian
RGCRFWidth_bipolarUnits = RGCRFWidth / bipolarSpacing;
RGCinputs2D = fspecial('gaussian', ceil(sqrt(nBipolars)), RGCRFWidth_bipolarUnits);
Tbp = TwoD2OneDTransform; %in units of bipolars
Tbp.rows = ceil(sqrt(nBipolars));
Tbp.cols = ceil(sqrt(nBipolars));
temp = Tbp.invert(RGCinputs2D);
RGC.inputGain = temp(1:nBipolars);

%set up model
M = FFModel;
M.layers = [bipolarLayer, RGCLayer];
M.nTrials = nTrials;

disp('done with model init');
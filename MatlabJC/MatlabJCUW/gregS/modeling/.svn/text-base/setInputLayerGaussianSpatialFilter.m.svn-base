function layer = setInputLayerGaussianSpatialFilter(layer,params)
units = layer.subunits;
%parameters
stimSize = params.stimSize; %square, (pixels)
bipolarSpacing = params.bipolarSpacing;
bipolarRFWidth = params.bipolarRFWidth;

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

disp(['Adding filters to ' num2str(nBipolars) ' bipolars...']);
for i=1:nBipolars;    
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
    units(i).filter = T.invert(f);
end
disp('done');

layers.subunits = units;


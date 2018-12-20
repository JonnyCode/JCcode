function L = makeConeSpatialSamplingLayer(params)
%parameters
stimSize = params.stimSize; %square, (pixels)
coneSpacing = params.coneSpacing;
coneRFWidth = params.coneRFWidth;

%set up cones
%Generate hexagonal grid
gridSize = ceil((stimSize./coneSpacing) * 1.1);
if rem(gridSize,2) == 1, gridSize = gridSize+1; end

Rad3Over2 = sqrt(3) / 2;
[X Y] = meshgrid(1:1:gridSize);
n = size(X,1);
X = Rad3Over2 * X;
Y = Y + repmat([0 0.5],[n,n/2]);

%set spacing
X = X * coneSpacing;
Y = Y * coneSpacing;

Ind = (X-coneRFWidth <= stimSize & Y-coneRFWidth <= stimSize);
X = X(Ind);
Y = Y(Ind);

%X,Y are now coordinates of cones centers
coneCenters = [X, Y];
nCones = numel(X);

%make gaussian 2D window
G = fspecial('gaussian',coneRFWidth*6,coneRFWidth);

disp(['Making ' num2str(nCones) ' cones...']);
for i=1:nCones;
    C(i) = Subunit;
    C(i).inputGain = [];
    %C(i).outputNL = @(x,T)sum(x); %make outputNL cone noise model
    
    %make filter
    f = zeros(stimSize,stimSize);
    startPointX = ceil(coneCenters(i,1)-coneRFWidth*3)+1; %3s.d of each side
    endPointX = floor(coneCenters(i,1)+coneRFWidth*3);    
    startPointY = ceil(coneCenters(i,2)-coneRFWidth*3)+1; %3s.d of each side
    endPointY = floor(coneCenters(i,2)+coneRFWidth*3);
    
    %edge effects
    startPoint_correctedX = max([startPointX,1]);
    endPoint_correctedX = min([endPointX,stimSize]);
    startPoint_correctedY = max([startPointY,1]);
    endPoint_correctedY = min([endPointY,stimSize]);
    
    fIndX = startPoint_correctedX:endPoint_correctedX;
    fIndY = startPoint_correctedY:endPoint_correctedY;
    GIndX = startPoint_correctedX - startPointX + 1:coneRFWidth*6-(endPointX - endPoint_correctedX);
    GIndY = startPoint_correctedY - startPointY + 1:coneRFWidth*6-(endPointY - endPoint_correctedY); 
    GIndX = GIndX(1:length(fIndX));
    GIndY = GIndY(1:length(fIndY));
    f(fIndX,fIndY) = G(GIndX,GIndY); 
    C(i).spatialPreFilter = f;    
end
disp('done');

%cone layer
L = Layer;
L.subunits = C;
%set transform
Ttime = FFT1DTransform;
L.transform = Ttime;




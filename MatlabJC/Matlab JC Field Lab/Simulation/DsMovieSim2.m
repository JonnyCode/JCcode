% script will run a set of simulated DS filters on a movie
% modified from 'DsMovieSim' to using more dynamic image motion and focus
% on issues in paper

% JC 2018/5/31

%% simulation

% load images
load('/Volumes/lab/Documents/Natural_Images/VhNaturalImages/VhImage_imk03002.mat')

% image/motion params
MovementType = 'jitter' ; % drift or jitter
SquareSize = true ; % zero=no square
SpeedStim = 8 ; % (pix/frame)
VhImage_size = size(VhImage) ;
FrameNumber = 5000 ; % how many frames to make the movie
DispLength = [400,400] ; % Y,X length of displayed image
VhImage_center = VhImage_size/2 ; %
MaskLength = VhImage_size-DispLength ;

if strcmp(MovementType,'drift') ;
    Delta = 50 ;
    Xnodes = [0:Delta:MaskLength(1)] ; % all possible stop/turn points
    Ynodes = [0:Delta:MaskLength(2)] ; 
    XnodesOrder = Xnodes(randi(length(Xnodes),[1,FrameNumber])) ; % list of nodes (only using a subset)
    YnodesOrder = Ynodes(randi(length(Ynodes),[1,FrameNumber])) ;
    DistPix = sqrt(diff(XnodesOrder).^2+diff(YnodesOrder).^2) ; % (pix) between stops
    DistFrames = DistPix/SpeedStim ;
    
    X=XnodesOrder(1) ;
    Y=YnodesOrder(1) ; % start at first node
    
    f=1 ; % set f
    while length(X)<=FrameNumber ;  % if all frames are not full
        Xtemp = interp1([1,DistFrames(f)],[XnodesOrder(f),XnodesOrder(f+1)],[1:DistFrames(f)]) ;
        Ytemp = interp1([1,DistFrames(f)],[YnodesOrder(f),YnodesOrder(f+1)],[1:DistFrames(f)]) ;
        X = [X,Xtemp] ;
        Y = [Y,Ytemp] ;
        f=f+1 ;
    end
    X = round(X(1:FrameNumber)) ;
    Y = round(Y(1:FrameNumber)) ;
    
    Mov = nan(DispLength(1),DispLength(2),FrameNumber) ;
    for f=1:FrameNumber ;   
        Mov(:,:,f) = VhImage([(X(f)+1):(X(f)+DispLength(1))],[(Y(f)+1):(Y(f)+DispLength(2))]) ; 
        %image(squeeze(Mov(f,:,:))); pause
    end
end



if strcmp(MovementType,'jitter') ;
    jitterStd = 60 ;
    jitterSmoothFrames = 300 ;
    MaxX = DispLength(1)/2 ; 
    MaxY = DispLength(2)/2 ; 
    
    X = normrnd(0,1,1,jitterStd) ; % drawl from norm distribution
    Y = normrnd(0,1,1,jitterStd) ;
    X = smooth(X,jitterSmoothFrames) ; % smooth
    Y = smooth(Y,jitterSmoothFrames) ;
    X(X<-MaxX)= -MaxX ; % rectify
    X(X>MaxX)= MaxX ;
    Y(Y<-MaxY)= -MaxY ; % rectify
    Y(Y>MaxY)= MaxY ;
    
    if SquareFlag ;
        Mov = zeros(DispLength(1),DispLength(2),FrameNumber) ;
        for f=1:FrameNumber ; 
            Xstrt = max(X(f) + MaxX - SquareSize,0) ;
            Xstrt = min(X(f) + MaxX,DispLength(1)) ;
            
            Ystrt = max(Y(f) + MaxY - SquareSize,0) ;
            Yend = min(Y(f) + MaxY,DispLength(2)) ;
            
            Mov(Xstrt:Xstrt,Ystrt:Yend,f) = 1 ;
            
            image(squeeze(Mov(:,:,f))); pause
        end
    end
end
        
Xdiff = [0,diff(X)] ;
Ydiff = [0,diff(Y)] ;
DirVect = cart2pol(Xdiff,Ydiff)*180/pi ;
DirVect(DirVect<0) = 360+DirVect(DirVect<0) ;

%% model cells and responses

% mosaic positions
DSdensity = 80 ; % /mm^2 (from rabbit Vaney 1994 ~20-50)
exclusion_mean = 150/4 ; % dendrite length (um)/(um/pix)
exclusion_sigma = 20/4 ; 
grid_size = [DispLength(1),DispLength(2)]-(exclusion_mean*2) ; % subtracting exclusion mean to keep cells of edges
num_cells = ceil((prod(grid_size)*4^2/1000^2)*DSdensity) ; % pix size of grid*(um/pix)/(um/mm)*density/mm^2

% directions and speeds
DsDir = [0:90:350] ; % perfered direction of cells (eg. [0:90:350,nan],[nan,nan,nan,nan,nan])
speed = SpeedStim ; % (pix/frame) speed preference of filters

% trf params
TrfTime = [1:ceil(ceil(exclusion_mean*2)/speed)] ; % this is the length of the ds filter
tpeak = [5,5,5,5] ;
peakRise = [100,100,100,100] ;
peakAmp = [1,1,1,1] ;
ttrough = [5,5,5,5] ;
troughDecay = [1,1,1,1] ;
troughAmp = [0,0,0,0] ;

% positions of DS cells
%    DsLocsSingle = make_serial_exclusion_mosaic(grid_size, num_cells, exclusion_mean, exclusion_sigma) ; %
%    DsLocsSingle = floor(DsLocsSingle)+ceil(exclusion_mean) ; % center

for DsType = 1:length(DsDir) ;
     DsLocs{DsType} = make_serial_exclusion_mosaic(grid_size, num_cells, exclusion_mean, exclusion_sigma) ; %
     DsLocs{DsType} = floor(DsLocs{DsType})+ceil(exclusion_mean) ; % center
%    DsLocs{DsType} = DsLocsSingle ;
    Params.centers = DsLocs{DsType} ; % (pix) X,Y positions of filter centers (Nx2), 
    Params.direction_preferences = ones(num_cells,1)*DsDir(DsType) ; % (deg) direction preference of each filter (Nx1)
    Params.speed = speed ; % (pix/frame) speed preference of filters
    Params.SrfFlag = true ; 
    Params.radi = exclusion_mean;
    %Params.Trf = simFilter(TrfTime,tpeak(DsType),peakRise(DsType),...
        %peakAmp(DsType),ttrough(DsType),troughDecay(DsType),troughAmp(DsType)) ;
    
    Lp{DsType} = LinDirFilters(Mov, Params) ;
end

% square Lp
for DsType = 1:size(Lp,2) ;
    Lp2{DsType} = Lp{DsType}.^2 ; % nonlinear
    %Lp2{DsType} = abs(Lp{DsType}) ;
end
    
% % normalize response for each image (contrast normalization)
% for im = 1:numImages; % for each image
%     imi = find(ImVect==imageNum(im)) ;
%     for DsType = 1:size(Lp,2) ;
%         for cells=1:num_cells ; % for each cell
%             Lp2{DsType}(cells,imi) = Lp2{DsType}(cells,imi)/median(Lp2{DsType}(cells,imi)) ;
%         end
%     end
% end

% rectify Lp
Thresh = 0.00 ; % fraction of max below which = 0 
Lp3 = Lp2 ;
for DsType = 1:size(Lp,2) ;
    for cells=1:num_cells ; % for each cell
        ri = Lp2{DsType}(cells,:)<max(Lp2{DsType}(cells,:))*Thresh ;
        Lp3{DsType}(cells,ri) = 0;
    end
end

%% correlation
for DsType = 1:size(Lp,2) ; % for each type
    for cells=1:num_cells ; % for each cell
        for DsType2 = 1:size(Lp,2) ; % for each type
            for cells2=1:num_cells ; % for each cell
                Corr{DsType}{DsType2}(cells,cells2) = corr(Lp3{DsType}(cells,:)',Lp3{DsType2}(cells2,:)') ;
                Dist{DsType}{DsType2}(cells,cells2) = pdist([DsLocs{DsType}(cells,:);DsLocs{DsType2}(cells2,:)]) ;
            end
        end
    end
end

%% OLE train and test
for TEMP_LOOP = 1:20 ;

% select cells
PopSubsetFlag = false ; % if false use all cells
NumberOleCells = 1 ; % number per type 
RandLocFlag = false ; % if true choose cells at random locations
RspSmoothFlag = false ; % smooth responses in ole
NumSmoothBins = 1 ; % number of frames to smooth over
SingleCellFlag = false ; % use only one cell to test OLE
AssumeWeightsFlag = true ; % assume wieghts by direction preference of cell

if PopSubsetFlag ;
    % cells by location
    if RandLocFlag ;
        for DsType = 1:length(DsDir) ;
            RandOrder = randperm(num_cells) ;
            OleCellsi{DsType} = RandOrder(1:NumberOleCells) ;
        end

    else % if it is local
        SeedSpot = [randi(DispLength(1)),randi(DispLength(2))] ; % [DispLength(1)/2,DispLength(2)/2] ;

        for DsType = 1:length(DsDir) ;
            for cells=1:num_cells ; % for each cell
                TempDist(cells) = sqrt(sum((SeedSpot-DsLocs{1}(cells,:)).^2)) ;
            end
            [~,DistOrder] = sort(TempDist) ;
            OleCellsi{DsType} = DistOrder(1:NumberOleCells) ; 
        end
    end
else
    for DsType = 1:length(DsDir) ;
        OleCellsi{DsType} = [1:num_cells] ; 
    end
end

% entire movie response mat
clear RspMat
c=1 ;
for DsType = 1:size(Lp3,2) ; % for each type
    for cells=1:length(OleCellsi{DsType}) ; % for each cell
        RspMat(:,c) = Lp3{DsType}(OleCellsi{DsType}(cells),:) ;
        c=c+1 ;
    end
end

if RspSmoothFlag ;
    for cells=1:size(RspMat,2) ; % for each cell
        RspMat(:,cells) = smooth(RspMat(:,cells),NumSmoothBins) ;
        ShiftPnts = floor(NumSmoothBins/2) ; % shift pnts so that smooth bin is causal
        RspMat(NumSmoothBins:end,cells) = RspMat(ShiftPnts+1:end-ShiftPnts,cells) ; % leaves first set innacurate
    end
end

ErrorHistX = [0:3:180] ;
clear Error

Traini = [1:floor(FrameNumber*3/4)] ; % ;
Testi =  [floor(FrameNumber*3/4)+1:FrameNumber];

RspMatTrain = RspMat(Traini,:) ;
DirVectTrain = DirVect(Traini) ;

% train OLE
clear OleWeights
if AssumeWeightsFlag ;
    c=1 ;
    for DsType = 1:size(Lp3,2) ; % for each type
        for cells=1:length(OleCellsi{DsType}) ; % for each cell
            [OleWeights(c,1),OleWeights(c,2)] = pol2cart(DsDir(DsType)*pi/180,1) ;
            c=c+1 ;
        end
    end
else
    [OleWeights, OlePolarWeights] = OleFinder(DirVectTrain,RspMatTrain) ;
end

% test OLE
RspMatTest = RspMat ;
DirVectTest = DirVect(Testi) ;

if SingleCellFlag ; % if you want only a single responsive cell
    for tp = 1:FrameNumber ; % for every point
        [mx,mxi] = max(RspMatTest(tp,:)) ; % find max
        RspMatTest(tp,:) = RspMatTest(tp,:)*0 ; %  zero all
        RspMatTest(tp,mxi) = mx ; % add only cell with max rate
    end
end
      
[dEst,dEstMag] = OleTester(OleWeights,RspMatTest) ;

% error analysis
for f=1:FrameNumber ; % for each estimated frame direction
    Error(f) = acuteAngle(DirVect(f),dEst(f)) ; % angle difference
end

ErrorMed = median(Error(Testi)) ; % median error
TEMP_LOOP_ErrorMed(TEMP_LOOP) = ErrorMed ;
    
end % END TEMP LOOP !!!!!!!!!!!

ErrorHist = hist(Error(Testi),ErrorHistX) ; % error histogram in test

% direction change triggered error
dirChangeThreshold = 15 ;
ChangeTrigErrorSize = 10 ; % number of Frames on each side of change

DirVectChangei = find(abs(diff(DirVect))>dirChangeThreshold) ; 

ChangeTrigError = zeros(1,ChangeTrigErrorSize*2+1) ;
Cnt = 0 ;
for tempi = 1:length(DirVectChangei)
    if DirVectChangei(tempi)>ChangeTrigErrorSize && ...
            DirVectChangei(tempi)<FrameNumber-ChangeTrigErrorSize ; % if you can get a full wave
        ChangeTrigError = ChangeTrigError + ...
            Error(DirVectChangei(tempi)-ChangeTrigErrorSize:DirVectChangei(tempi)+ChangeTrigErrorSize) ;
        Cnt = Cnt+1 ;
    end
end
ChangeTrigError = ChangeTrigError/Cnt ; % average

% RspMat stats    
RspMat_sum = sum(RspMat') ; % sum all neurons

RspMat_FracAbove0 = sum(RspMat_sum>0)/length(RspMat_sum) ;
RspMat_NumCells = sum(RspMat'>0) ;

% median error per number responsive cells
RspMat_NumCellsX = unique(RspMat_NumCells) ;
for nc = 1:length(RspMat_NumCellsX) ; % for each number 
    Error_FunNumCells(nc) = median(Error(RspMat_NumCells==RspMat_NumCellsX(nc))) ; % median error in bins with X# cells
end

%% figures
Color_list = {'r','c','k','b'} ;

figure % correlation
for DsType = 1:4 ; % for each type
    for cells=1:size(Lp{DsType},1) ; % for each cell
        for DsType2 = 1:4 ; % for each type
            for cells2=1:size(Lp{DsType2},1) ; % for each cell
                if DsType == DsType2 ; 
                    PlotColor = 'ko' ;
                elseif abs(DsType-DsType2)==1 | abs(DsType-DsType2)==3 ;
                    PlotColor = 'ro' ;
                else
                    PlotColor = 'bo' ;
                end
                
                plot(Dist{DsType}{DsType2}(cells,cells2),Corr{DsType}{DsType2}(cells,cells2),PlotColor)
                hold on
            end
        end
    end
end
xlabel('dist (pix)')
ylabel('corr coef')
    
figure % OLE estimate
plot(DirVect,'k')
hold on
plot(dEst,'r')
plot(Testi,dEst(Testi),'g')
xlabel('frame')
ylabel('direction (deg)')

figure ; % change triggered error
plot(ChangeTrigError,'g-')
xlabel('frames')
ylabel('error (deg)')

figure ; % plot mosaic
for DsType = 1:length(DsDir) ;
    for cells = 1:num_cells ;
        plot(DsLocs{DsType}(cells,1),DsLocs{DsType}(cells,2),[Color_list{DsType},'o'])
        hold on
        %drawCircle(DsLocs{DsType}(cells,1),DsLocs{DsType}(cells,2),exclusion_mean,'LineStyle','-','color',Color_list{DsType})
    end
end
for DsType = 1:length(DsDir) ;
    for cells = 1:length(OleCellsi{DsType}) ;
        plot(DsLocs{DsType}(OleCellsi{DsType}(cells),1),...
            DsLocs{DsType}(OleCellsi{DsType}(cells),2),'or')
    end
end
   
figure % errror hist
plot(ErrorHistX,ErrorHist)

% population figure from excel
% controlNumCells = [1,2,5,10,30,34] ;
% control = [] ;
% randomGrowth = [] ;
% singleCellGrowth = [] ;
% sameMosaic = [] ;
% control_std = [] ;
% randomGrowth_std = [] ;
% singleCellGrowth_std = [] ;
% sameMosaic_std = [] ;
% 
% denseNumCells = [1,2,5,10,30,60,136] ;
% denseMosaic = [] ;
% denseMosaic_std = [] ;

figure
errorbar(controlNumCells,control,control_std/sqrt(20))
hold on
errorbar(controlNumCells,randomGrowth,randomGrowth_std/sqrt(20))
errorbar(controlNumCells,singleCellGrowth,singleCellGrowth_std/sqrt(20))
errorbar(controlNumCells,sameMosaic,sameMosaic_std/sqrt(20))
errorbar(denseNumCells,denseMosaic,denseMosaic_std/sqrt(20))
xlabel('number of quads')
ylabel('median error')
hgsave('errorCellNumber')
print(gcf,'-dpdf','errorCellNumber')


%% NOT USED
RspCorrFlag = false ; % run code to check calculate correlation in rsp mat
if RspCorrFlag ; % correlation in RspMat -averaged across direction
    clear RspMati
    for DsType = 1:size(Lp3,2) ; % for each type
        RspMati(DsType,:) = [1:NumberOleCells]+(DsType-1)*NumberOleCells ; % index for each type within RspMat
    end

    CorrMat = zeros(NumberOleCells*4) ; % prep mat
    for tempi = 2:length(DirVectChangei) ; % for each direction (skip first and last for convenience)
        CorrMat = CorrMat+corrcoef(RspMat([DirVectChangei(tempi-1):DirVectChangei(tempi)],:)) ;
    end
    CorrMat = CorrMat/(length(DirVectChangei)-1) ;

    RspMatCorr = zeros(size(Lp,2)) ;
    RspMatCorrDenom = zeros(size(Lp,2)) ;

    for DsType = 1:size(Lp,2) ; % for each type
        for DsType2 = 1:size(Lp,2) ; % for each type
            for cells=1:NumberOleCells ; % for each cell
                for cells2=1:NumberOleCells ; % for each cell
                    if ~(DsType==DsType2 && cells==cells2) ; % if its not the same cell
                        RspMatCorr(DsType,DsType2) = RspMatCorr(DsType,DsType2) ...
                            + CorrMat(RspMati(DsType,cells),RspMati(DsType2,cells2)) ; % sum correlation
                        RspMatCorrDenom(DsType,DsType2) = RspMatCorrDenom(DsType,DsType2)+1 ; % count
                    end
                end
            end
        end
    end
    RspMatCorr =  RspMatCorr./RspMatCorrDenom ;               

    RspMatCorrIntra(TEMP_LOOP) = mean(diag(RspMatCorr)) ; % average correlation within type      
    upt = (triu(ones(4),1)) ; % average correlation across types
    RspMatCorrInter(TEMP_LOOP) = mean(RspMatCorr(triu(ones(4),1)==1)) ;
end



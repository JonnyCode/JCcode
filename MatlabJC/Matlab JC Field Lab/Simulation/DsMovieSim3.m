% script will run a set of simulated DS filters on a movie
% modified from 'DsMovieSim2' to use jitter square or image
% changed how 'jitter' calculates vectors 
% added 'steps' (like 'drift' but calculated differently)

% JC 2018/5/31

%% simulation

% load images
%load('/Volumes/lab/Documents/Natural_Images/VhNaturalImages/VhImage_imk03002.mat')
load('/Volumes/lab/Documents/Natural_Images/VhNaturalImages/VhImage_imk00002.mat')

% image/motion params
MovementType = 'steps' ; % drift or jitter
SquareSize = 0 ; % zero=no square
SpeedStim = 8 ; % (pix/frame)
FrameNumber = 2000 ; % how many frames to make the movie
DispLength = [600,600] ; % Y,X length of displayed image
rng(1) ; % sets the seed for random number generator so if framenumber is the same than trajectory will be the same
ImageSizeX = size(VhImage,1) ; % hieght
ImageSizeY = size(VhImage,2) ; % width

if SquareSize>0;
    MaxX = DispLength(1)-SquareSize ;
    MaxY = DispLength(2)-SquareSize ;
else
    MaxX = ImageSizeX - DispLength(1) ;
    MaxY = ImageSizeY - DispLength(2) ;
end
    

if strcmp(MovementType,'jitter') ;
    jitterStd = 8 ;
    
    X(1) = DispLength(1)/2 ; % initial position 
    Y(1) = DispLength(2)/2 ;
    
    for f=2:FrameNumber ;
        X(f) = X(f-1)+round(normrnd(0,jitterStd,1,1)) ; % take step
        
        if X(f)>MaxX ;
            X(f) = DispLength(1)-SquareSize ;
        end
        if X(f)<1 ;
            X(f) = 1 ;
        end
        
        Y(f) = Y(f-1)+round(normrnd(0,jitterStd,1,1)) ; % take step
        
        if Y(f)>MaxY ;
            Y(f) = DispLength(2)-SquareSize ;
        end
        if Y(f)<1 ;
            Y(f) = 1 ;
        end
    end
end

if strcmp(MovementType,'steps') ;
    
    X(1) = DispLength(1)/2 ; % initial position 
    Y(1) = DispLength(2)/2 ; 
    xdiff = randi([-SpeedStim SpeedStim],1) ; % initial direction
    ydiff = randi([-SpeedStim SpeedStim],1) ;
    
    for f=2:FrameNumber ; % for each frame
        while (X(f-1)+xdiff)>MaxX || (X(f-1)+xdiff)<1 || (ydiff+xdiff)==0 ...
               || (Y(f-1)+ydiff)>MaxY || (Y(f-1)+ydiff)<1 % hit edge or not moving
                
            xdiff = randi([-SpeedStim SpeedStim],1) ; % pick new directions
            ydiff = randi([-SpeedStim SpeedStim],1) ;
        end
        
        X(f) = X(f-1)+xdiff ; % go in that direction
        Y(f) = Y(f-1)+ydiff ;
    end
end


% make movie
Mov = zeros(DispLength(1),DispLength(2),FrameNumber) ; % prep mat

if SquareSize>0;
    for f=1:FrameNumber ; 
        Xstrt = X(f) ;
        Xend = X(f)+SquareSize-1 ;

        Ystrt = Y(f) ;
        Yend = Y(f)+SquareSize-1 ;

        Mov(Xstrt:Xend,Ystrt:Yend,f) = 1 ;

        %imagesc(squeeze(Mov(:,:,f))); pause
    end
else % if its a image
    for f=1:FrameNumber ; 
        Mov(:,:,f) = VhImage([(X(f)+1):(X(f)+DispLength(1))],[(Y(f)+1):(Y(f)+DispLength(2))]) ; 
        %imagesc(squeeze(Mov(:,:,f))); pause
    end
end 


% direction of motion
Xdiff = [0,diff(X)] ;
Ydiff = [0,diff(Y)] ;
DirVect = cart2pol(Xdiff,Ydiff)*180/pi ;
DirVect(DirVect<0) = 360+DirVect(DirVect<0) ;


%% model cells and responses

% mosaic positions
DSdensity = 30 ; % /mm^2 (from rabbit Vaney 1994 ~20-50)
exclusion_mean = 150/4 ; % dendrite length (um)/(um/pix)
exclusion_sigma = 20/4 ; 
grid_size = [DispLength(1),DispLength(2)]-(exclusion_mean*2) ; % subtracting exclusion mean to keep cells of edges
num_cells = ceil((prod(grid_size)*4^2/1000^2)*DSdensity) ; % pix size of grid*(um/pix)/(um/mm)*density/mm^2

% directions and speeds
DsDir = [0:90:350] ; % perfered direction of cells (eg. [0:90:350,nan],[nan,nan,nan,nan,nan])
speed = SpeedStim ; % (pix/frame) speed preference of filters

% trf params
TrfTime = [1:ceil(ceil(exclusion_mean*2)/speed)] ; % this is the length of the ds filter
tpeak = [4,4,4,4] ;
peakRise = [100,100,100,100] ;
peakAmp = [1,1,1,1] ;
ttrough = [5,5,5,5] ;
troughDecay = [1,1,1,1] ;
troughAmp = [0,0,0,0] ;

for DsType = 1:length(DsDir) ;
     DsLocs{DsType} = make_serial_exclusion_mosaic(grid_size, num_cells, exclusion_mean, exclusion_sigma) ; %
     DsLocs{DsType} = floor(DsLocs{DsType})+ceil(exclusion_mean) ; % center
end

for DsType = 1:length(DsDir) ;
    Params.centers = DsLocs{DsType} ; % (pix) X,Y positions of filter centers (Nx2), 
    Params.direction_preferences = ones(num_cells,1)*DsDir(DsType) ; % (deg) direction preference of each filter (Nx1)
    Params.speed = speed ; % (pix/frame) speed preference of filters
    Params.SrfFlag = true ; 
    Params.radi = exclusion_mean;
    %Params.Trf = simFilter(TrfTime,tpeak(DsType),peakRise(DsType),...
        %peakAmp(DsType),ttrough(DsType),troughDecay(DsType),troughAmp(DsType)) ;
    
    Lp{DsType} = LinDirFilters(Mov, Params) ;
end

% on-off Lp
for DsType = 1:size(Lp,2) ;
    Lp2{DsType} = abs(Lp{DsType}) ;
end

% rectify Lp
Thresh = 0.085 ; % fraction of max below which = 0 (0.275 provides similar response sparsity to data)
TrueThresh = max(Lp2{DsType}(:))*Thresh ; % threshold for all cells
for DsType = 1:size(Lp,2) ;
    for cells=1:num_cells ; % for each cell
        Lp3{DsType}(cells,:) = Lp2{DsType}(cells,:) - TrueThresh ; % subtract
        ri = Lp3{DsType}(cells,:)<0; % zero below threshold
        Lp3{DsType}(cells,ri) = 0;
    end
end

%% respones sparity measure
AllRsp = cell2mat(Lp3) ; % all ds cells
FractionNoResponse = sum(AllRsp(:)==0)/length(AllRsp(:))

%% save responses and movie

SavePath = '/Users/jcafaro/Documents/AnalysisFigures/NatStimDs/2018Paper/DsMovieSim3/mfile/' ;
RunId = datestr(now,'yyyymmddhhMM') ;

save([SavePath,'Mov2SimOnly_RunId',num2str(RunId)],'-v7.3') ; %save all


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
NumOleCellsSet = [2.^[0:7],num_cells]; % number of cells to test

for NumCells_loop = 1:length(NumOleCellsSet) ; % for each number of cells in decoder

    % cells select params
    PopSubsetFlag = true ; % if false use all cells
    NumberOleCells = NumOleCellsSet(NumCells_loop) ; % number per type 
    RandLocFlag = false; % if true choose cells at random locations
    RspSmoothFlag = false ; % smooth responses in ole
    NumSmoothBins = 1 ; % number of frames to smooth over
    SingleCellFlag = false ; % use only one cell to test OLE
    AssumeWeightsFlag = true ; % assume wieghts by direction preference of cell
    NumReps = 50 ; % number of repeats
    
    for Rep_loop = 1:NumReps ; % repeat with different subsets of cells

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
                        TempDist(cells) = sqrt(sum((SeedSpot-DsLocs{DsType}(cells,:)).^2)) ;
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
        ErrorMed_RepLoop(Rep_loop) = ErrorMed ;

    end % END Rep LOOP

ErrorMed_NumCellsLoop(NumCells_loop) = mean(ErrorMed_RepLoop) ; % average across repeats
ErrorMed_NumCellsLoop_Sem(NumCells_loop) = std(ErrorMed_RepLoop)/sqrt(NumReps) ; % std across repeats
end


%% save
SavePath = '/Users/jcafaro/Documents/AnalysisFigures/NatStimDs/2018Paper/DsMovieSim3/mfile/' ;
SaveFigPath = '/Users/jcafaro/Documents/AnalysisFigures/NatStimDs/2018Paper/DsMovieSim3/figs/' ;
RunId = datestr(now,'yyyymmddhhMM') 

clear Mov % save room
save([SavePath,'RunId',num2str(RunId)],'-v7.3') ; %save all

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
%plot(Testi,dEst(Testi),'g')
xlabel('frame')
ylabel('direction (deg)')

saveas(gcf,[SaveFigPath,'DsMovieSim3_RunId',num2str(RunId)])
print(gcf, '-dpdf',[SaveFigPath,'DsMovieSim3_RunId',num2str(RunId)])

figure ; % change triggered error
plot(ChangeTrigError,'g-')
xlabel('frames')
ylabel('error (deg)')

figure ; % plot full mosaic
for DsType = 1:length(DsDir) ;
    for cells = 1:num_cells ;
        plot(DsLocs{DsType}(cells,1),DsLocs{DsType}(cells,2),[Color_list{DsType},'o'])
        hold on
        drawCircle(DsLocs{DsType}(cells,1),DsLocs{DsType}(cells,2),exclusion_mean,'LineStyle','-','color',Color_list{DsType})
    end
end
for DsType = 1:length(DsDir) ;
    for cells = 1:length(OleCellsi{DsType}) ;
        plot(DsLocs{DsType}(OleCellsi{DsType}(cells),1),...
            DsLocs{DsType}(OleCellsi{DsType}(cells),2),'or')
    end
end
   

figure ; % plot mosaic of cells decoded
for DsType = 1:length(DsDir) ;
    for cells = 1:length(OleCellsi{DsType}) ; % for each cell decoded
        %plot(DsLocs{DsType}(OleCellsi{DsType}(cells),1),DsLocs{DsType}(OleCellsi{DsType}(cells),2),[Color_list{DsType},'o'])
        hold on
        drawCircle(DsLocs{DsType}(OleCellsi{DsType}(cells),1),DsLocs{DsType}(OleCellsi{DsType}(cells),2),exclusion_mean,'LineStyle','-','color',Color_list{DsType})
    end
end
for DsType = 1:length(DsDir) ;
    for cells = 1:length(OleCellsi{DsType}) ;
        plot(DsLocs{DsType}(OleCellsi{DsType}(cells),1),...
            DsLocs{DsType}(OleCellsi{DsType}(cells),2),'or')
    end
end

axis([0 600 0 600])
saveas(gcf,[SaveFigPath,'MosaicLocal_RunId',num2str(RunId)])
print(gcf, '-dpdf',[SaveFigPath,'MosaicLocal_RunId',num2str(RunId)])


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


figure ; % movie of response
for f=1:FrameNumber ;
    for DsType = 1:length(DsDir) ;
        for cells = 1:num_cells ;
            plot(DsLocs{DsType}(cells,1),DsLocs{DsType}(cells,2),[Color_list{DsType},'o'],'MarkerSize',.1+10*Lp3{DsType}(cells,f)/max(Lp3{DsType}(:)))
            hold on
            %drawCircle(DsLocs{DsType}(cells,1),DsLocs{DsType}(cells,2),exclusion_mean,'LineStyle','-','color',Color_list{DsType})
        end
    end
    hold off
    pause
end

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

    RspMatCorrIntra(Rep_loop) = mean(diag(RspMatCorr)) ; % average correlation within type      
    upt = (triu(ones(4),1)) ; % average correlation across types
    RspMatCorrInter(Rep_loop) = mean(RspMatCorr(triu(ones(4),1)==1)) ;
end



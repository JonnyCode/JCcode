% script will run a set of simulated DS filters on a movie

% JC 2017/10/17

%% simulation

imageNum = [1:13] ; % image numbers to test
numImages = length(imageNum) ;

% load image slide movie
load('/Users/jcafaro/Documents/MATLAB/JCcode/MatlabJC/Matlab JC Field Lab/Simulation/DsMoviesForSim/Image_Slide_Stim2.mat')
figure ; imagesc(FullFrame{imageNum(1)}{1,1}) ; % fullFrame{image}{direction,frame}

Movie_dir = repmat([0:45:350],[numImages,1]) ; 
FullFrameXRange = [300:500] ; % selected x points (to exclude mask)
FullFrameYRange = [200:400] ;

numd = size(FullFrame{imageNum(1)},1) ; % number of directions
numf = size(FullFrame{imageNum(1)},2) ; % number of frames

% put into 1 image into Movie mat
Mov = nans(length(FullFrameYRange),length(FullFrameXRange),numf*numImages) ;
t=1 ;
for im = 1:length(imageNum) ;  % for each image 
    for d = 1:numd ; % for each direction
        for f = 1:numf ; % for each frame
            Mov(:,:,t) = FullFrame{imageNum(im)}{d,f}(FullFrameYRange,FullFrameXRange) ;
            t=t+1 ;
        end
    end
end

clear FullFrame % clear for memmory

% entire dir mat
DirVect = [] ;
ImVect = [] ;
for im = 1:numImages; % for each image
    for d = 1:size(Movie_dir,2) ; % for each direction
        DirVect = [DirVect;repmat(Movie_dir(im,d),[numf,1])] ;
        ImVect = [ImVect; repmat(imageNum(im),[numf,1])] ;
    end
end

% mosaic positions
DSdensity = 30 ; % /mm^2 (from rabbit Vaney 1994 ~20-50)
exclusion_mean = 150/4 ; % dendrite length (um)/(um/pix)
exclusion_sigma = 20/4 ; 
grid_size = [size(Mov,1),size(Mov,2)]-(exclusion_mean*2) ; % subtracting exclusion mean to keep cells of edges
num_cells = ceil((prod(grid_size)*4^2/1000^2)*DSdensity) ; % pix size of grid*(um/pix)/(um/mm)*density/mm^2

% directions and speeds
DsDir = [0:90:350] ; % perfered direction of cells (eg. [0:90:350,nan],[nan,nan,nan,nan,nan])
speed = 8 ; % (pix/frame) speed preference of filters

% trf params
TrfTime = [1:ceil(ceil(exclusion_mean*2)/speed)] ; % this is the length of the ds filter
tpeak = [5,5,5,5,2,4,6,8] ;
peakRise = [100,100,100,100,3,3,3,3] ;
peakAmp = [1,1,1,1,1,1,1,1] ;
ttrough = [5,5,5,5,4,5,6,7] ;
troughDecay = [1,1,1,1,2,2,2,2] ;
troughAmp = [0,0,0,0,0,0,0,0] ;

% positions of DS cells
for DsType = 1:length(DsDir) ;
    DsLocs{DsType} = make_serial_exclusion_mosaic(grid_size, num_cells, exclusion_mean, exclusion_sigma) ; %
    DsLocs{DsType} = floor(DsLocs{DsType})+ceil(exclusion_mean) ; % center
    
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
    Lp2{DsType} = Lp{DsType}.^2 ;
end
    
% normalize response for each image (contrast normalization)
for im = 1:numImages; % for each image
    imi = find(ImVect==imageNum(im)) ;
    for DsType = 1:size(Lp,2) ;
        for cells=1:num_cells ; % for each cell
            Lp2{DsType}(cells,imi) = Lp2{DsType}(cells,imi)/median(Lp2{DsType}(cells,imi)) ;
        end
    end
end

% rectify Lp
Thresh = 0.2 ; % fraction of max below which = 0 
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
                Temp = xcov(Lp3{DsType}(cells,:),Lp3{DsType2}(cells2,:),'coef') ;
                corr{DsType}{DsType2}(cells,cells2) = Temp(ceil(length(Temp)/2)) ;
                Dist{DsType}{DsType2}(cells,cells2) = pdist([DsLocs{DsType}(cells,:);DsLocs{DsType2}(cells2,:)]) ;
            end
        end
    end
end
    
AllCov = corrcoef(cell2mat(Lp3')') ; % corr coef of all cells

%% OLE train and test
OqeFlag = false ;
ModifierFlag = false ; % add modifier to OLE

KillNonDsFlag = false ;
KillDsFlag = false ;

KillDsSingleTermFlag = false ;
KillNonDsSingleTermFlag = false ;
KillNonDsNonDsCrossTermFlag = false ;
KillDsNonDsCrossTermFlag = false ;
KillDsDsCrossTermFlag = false ;

% crossterm indicies
NumOleCells = num_cells*size(Lp,2) ; % number of cells

CrossTermCelli = [] ;
for cells = 2:NumOleCells ; 
    CrossTermCelli = [CrossTermCelli,[repmat(cells,[1,cells-1]);[1:cells-1]]] ; % cellA;cellB in order of crossterms
end
CrossTermCellTypes = ceil(CrossTermCelli/num_cells) ; % types of cells in each cross terms

% first term indicies
TermCellTypes = ceil([1:NumOleCells]/num_cells) ;

% set flagged terms to zero
NonDsTypesi = find(isnan(DsDir)) ; % types that are non-DS

CtKilli = [] ;
if KillNonDsSingleTermFlag ;
    CtKilli = [CtKilli,find(ismember(TermCellTypes,NonDsTypesi))] ;
end

if KillDsSingleTermFlag ;
    CtKilli = [CtKilli,find(~ismember(TermCellTypes,NonDsTypesi))] ;
end

if KillNonDsNonDsCrossTermFlag ;
    CtKilli = [CtKilli,find(sum(ismember(CrossTermCellTypes,NonDsTypesi),1)==2)+NumOleCells];
end

if KillDsNonDsCrossTermFlag ;
    CtKilli = [CtKilli,find(sum(ismember(CrossTermCellTypes,NonDsTypesi),1)==1)+NumOleCells] ;
end
 
if KillDsDsCrossTermFlag ;
    CtKilli = [CtKilli,find(sum(ismember(CrossTermCellTypes,NonDsTypesi),1)==0)+NumOleCells] ;
end

% entire movie response mat
clear RspMat
c=1 ;
for DsType = 1:size(Lp3,2) ; % for each type
    for cells=1:size(Lp3{DsType},1) ; % for each cell
        RspMat(:,c) = Lp3{DsType}(cells,:) ;
        c=c+1 ;
    end
end

ErrorHistX = [0:3:180] ;
clear Error
% hold out 1 image for test, train on other
for im = 1:length(imageNum) ;  % for each image to test 
    Traini = find(ImVect~=imageNum(im))  ; % [1:length(ImVect)] ;
    Testi = find(ImVect==imageNum(im)) ;
    
    RspMatTrain = RspMat(Traini,:) ;
    DirVectTrain = DirVect(Traini) ;
    
    if OqeFlag ;
        [OleWeights, OlePolarWeights] = OleFinder(DirVectTrain,RspMatTrain,'OQE') ;
    else
        [OleWeights, OlePolarWeights] = OleFinder(DirVectTrain,RspMatTrain) ;
    end

    % test OLE
    RspMatTest = RspMat(Testi,:) ;
    DirVectTest(im,:) = DirVect(Testi) ;
    
    c=1 ;
    for DsType = 1:size(Lp,2) ; % for each type
        for cells=1:size(Lp{DsType},1) ; % for each cell
            if KillNonDsFlag ; %
                if isnan(DsDir(DsType)) ;
                    RspMatTest(:,c) = RspMatTest(:,c)*0;
                end
            end
            if KillDsFlag ; % 
                if ~isnan(DsDir(DsType)) ;
                    RspMatTest(:,c) = RspMatTest(:,c)*0;
                end 
            end
            c=c+1 ;
        end
    end

    % modify wieghts if flagged
    OleWeightsTest = OleWeights ;
    if ~isempty(CtKilli) ;
        OleWeightsTest(CtKilli,:) = 0 ;
    end
    
    if OqeFlag ;
        [dEst(im,:),dEstMag(im,:)] = OleTester(OleWeightsTest,RspMatTest,'OQE') ;
    else
        [dEst(im,:),dEstMag(im,:)] = OleTester(OleWeightsTest,RspMatTest) ;
    end

    % error analysis
    for f=1:length(dEst) ; % for each estimated frame direction
        Error(im,f) = acuteAngle(DirVectTest(im,f),dEst(im,f)) ; % angle difference
    end
    
    ErrorHistbyIm(im,:) = hist(Error(im,:),ErrorHistX) ;
    
    if ModifierFlag ; % if you want to test the modifier 
        % alpha modifier (negative from other cells)
        alphaFact = 0.5 ;
        alphaMod = -alphaFact*(AllCov*RspMatTest.^2')'*OleWeightsTest(1:end-1,:) ; % modifier in X,Y
        
        % beta modifier (positive from product of cell and other cells)
        betaFact = 0.5 ;
        RWmat{1} = RspMatTest'.*repmat(OleWeightsTest(1:end-1,1),1,size(RspMatTest,1)) ;
        RWmat{2} = RspMatTest'.*repmat(OleWeightsTest(1:end-1,2),1,size(RspMatTest,1)) ;
        betaMod(:,1) = betaFact*diag((AllCov*RspMatTest')'*RWmat{1}) ;
        betaMod(:,2) = betaFact*diag((AllCov*RspMatTest')'*RWmat{2}) ;
        
        CartMod = alphaMod+betaMod ; % cartisian modifiers
        
        % polar estimate
        PolarMod = atan2d(CartMod(:,2), CartMod(:,1)) ; % four quadrant inverse tangent
        PolarMod(PolarMod<0) = PolarMod(PolarMod<0)+360 ; % no negatives
        PolarModMag = sqrt(CartMod(:,2).^2 + CartMod(:,1).^2) ;
        
        for imtb = 1:length(dEst(im,:)) ; % for each image time bin
            Temp = PolarVectorAddition([dEst(im,imtb),dEstMag(im,imtb); PolarMod(imtb),PolarModMag(imtb)]) ; 
            dEstMod(im,imtb) = Temp(1) ;
            dEstModMag(im,imtb) = Temp(2) ;
        end
        
        % error analysis
        for f=1:length(dEstMod) ; % for each estimated frame direction
            ErrorMod(im,f) = acuteAngle(DirVectTest(im,f),dEstMod(im,f)) ; % angle difference
        end

        ErrorModHistbyIm(im,:) = hist(ErrorMod(im,:),ErrorHistX) ;
    end
    
end

ErrorHist = hist(Error(:),ErrorHistX) ;

if ModifierFlag ; % if you want to test the modifier 
    ErrorHistMod = hist(ErrorMod(:),ErrorHistX) ;
    %ErrorHistMod_median = 
end

% mag hist
MagHist = hist(dEstMag(:),dEstMagHistX) ;

dEstMagHistX = [0:max(dEstMag(:))/100:max(dEstMag(:))] ; 
for b = 1:length(dEstMagHistX) ;  
    ErrorMedian(b) = median(Error(dEstMag>dEstMagHistX(b))) ;
end

%% crossterm analysis on OQE

% direction of cross terms 
for cells = 1:NumOleCells ; % for each cell in the OQE
    [temp,cti] = find(CrossTermCelli==cells) ; % find cross terms with that cell 
    cti = cti + NumOleCells ; % adjust for first set of terms which are not cross terms - index of OleWeights

    CrsTermVecSum(cells,:) = PolarVectorAddition(OlePolarWeights(cti,:)) ; % vector sum of cross terms

    TermCrsTermAngDiff(cells) = acuteAngle(CrsTermVecSum(cells,1),OlePolarWeights(cells,1)) ;
end
    
% angle and distance between cell locations
for cti = 1:size(CrossTermCelli,2) ; % for each cross term
    DsType1 = ceil(CrossTermCelli(1,cti)/size(Lp3{1},1)) ; % cell type 
    cell1 = CrossTermCelli(1,cti)-(DsType1-1)*size(Lp3{1},1) ; % cell #
    
    DsType2 = ceil(CrossTermCelli(2,cti)/size(Lp3{1},1)) ; % cell type 
    cell2 = CrossTermCelli(2,cti)-(DsType2-1)*size(Lp3{1},1) ; % cell #
    
    PosDelta = (DsLocs{DsType1}(cell1,:)-DsLocs{DsType2}(cell2,:)) ; % position delta
    CrossTermDist(cti) = sqrt(sum(PosDelta.^2)) ;
    CrossTermAngle(cti) = atan2d(PosDelta(1),PosDelta(2)) ;
    if CrossTermAngle(cti)<0 ; 
        CrossTermAngle(cti) = 360+CrossTermAngle(cti) ;
    end
end
    
%% figures

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
                
                plot(Dist{DsType}{DsType2}(cells,cells2),corr{DsType}{DsType2}(cells,cells2),PlotColor)
                hold on
            end
        end
    end
end
xlabel('dist (pix)')
ylabel('corr coef')
    
figure % OLE estimate
for im = 1:length(imageNum) ;  % for each image to test 
    plot(dEst(im,:))
    hold on
    plot(dEstMod(im,:),'g') % with modifier
    plot(DirVectTest(im,:),'r')
    pause
    hold off
end
xlabel('frame')
ylabel('direction (deg)')

figure ; %error
subplot(3,1,1)
for im = 1:length(imageNum) ;  % for each image to test 
    %plot(ErrorHistX,cumsum(ErrorHistbyIm(im,:)/sum(ErrorHistbyIm(im,:))),'c:')  
    hold on
end
plot(ErrorHistX,cumsum(ErrorHist/sum(ErrorHist)))
xlabel('error (deg)')
ylabel('fraction observed')

subplot(3,1,2)
plot(dEstMagHistX,cumsum(MagHist)/sum(MagHist))
hold on
xlabel('est mag')
ylabel('fraction observed')

subplot(3,1,3)
plot(dEstMagHistX,ErrorMedian)
hold on
xlabel('est mag')
ylabel('error median')



figure % cross terms wieghts - all cells
for cells = 1:NumOleCells ;
    [temp,cti] = find(CrossTermCelli==cells) ; % find cross terms with that cell
    cti = cti + NumOleCells ; % adjust for first set of terms which are not cross terms
    
    ctVs = PolarVectorAddition(OlePolarWeights(cti,:)) ; % vector sum of cross terms
    
    clf ;
    polar(0,max(OlePolarWeights(cti,2)),'b') ; % adjust axis so longest vector can be seen
    hold on
    for cells2 = 1:length(cti) ; % plot cross terms
        polar([0,OlePolarWeights(cti(cells2),1)*pi/180],[0,OlePolarWeights(cti(cells2),2)],'b')
    end
    polar([0,OlePolarWeights(cells,1)*pi/180],[0,OlePolarWeights(cells,2)],'r') ; % plot cell
    polar([0,ctVs(1)*pi/180],[0,ctVs(2)],'g') ; % plot cell
    pause
end

figure % cross terms vs single cell term histograms
hist(TermCrsTermAngDiff,[0:10:200])
xlabel('Angle Cross Term and single cell') 
ylabel('number of cells obs')

figure % cross terms wieghts - all cells
for cells = 1:NumOleCells ; % for each cell
    [temp,cti] = find(CrossTermCelli==cells) ; % find cross terms with that cell
    
    %plot(CrossTermAngle(cti),OlePolarWeights(cti+NumOleCells,1),'*')
    plot(CrossTermAngle(cti),OlePolarWeights(cti+NumOleCells,2),'*')
    hold on
    pause
end
 
figure % cross terms wieghts - all cells
for cti = 1:size(CrossTermCelli,2) ; % for each cross term

    DsType1 = ceil(CrossTermCelli(1,cti)/size(Lp3{1},1)) ; % cell type 
    cell1 = CrossTermCelli(1,cti)-(DsType1-1)*size(Lp3{1},1) ; % cell #

    DsType2 = ceil(CrossTermCelli(2,cti)/size(Lp3{1},1)) ; % cell type 
    cell2 = CrossTermCelli(2,cti)-(DsType2-1)*size(Lp3{1},1) ; % cell #
    
    subplot(1,2,1)
    polar([0,OlePolarWeights(CrossTermCelli(1,cti),1)*pi/180],...
        [0,1],'r') ; % cell 1
    hold on
    polar([0,OlePolarWeights(CrossTermCelli(2,cti),1)*pi/180],...
        [0,1],'b') ; % cell 2
    polar([0,OlePolarWeights(cti+NumOleCells,1)*pi/180],...
        [0,1],'k') ; % cross term
    hold off
    
    subplot(1,2,2)
    plot(DsLocs{DsType1}(cell1,1),DsLocs{DsType1}(cell1,2),'r*')
    hold on
    plot(DsLocs{DsType2}(cell2,1),DsLocs{DsType2}(cell2,2),'b*')
    axis([0 300 0 300])
    hold off
    legend(num2str(DsType1),num2str(DsType2))
    pause

end

    
function ForIgor = ImSliderConcatAnalysisV6(DataBlock, DB, Params)


% this function is adapted from 'ImSliderConcatAnalysisV5' to focus on
% multivariate variance and covariance 

% JC 3/23/2018 

LoadOldMatFileFlag = false; % true skip calculations and just load
PlotFigsFlag = true ;
saveFigPath = ['/Users/jcafaro/Documents/AnalysisFigures/'] ;
numElectrodeLayers = 2 ; % number of electrode layers surrounding the max that the center is calculated over

StimOnOffset = 1 ; % (sec) time after stim begins to start counting spikes
DgOffTime = 8 ; % (sec) time ds typing drifiting grating turned off
StimOffTime = .5 ; % (sec) time before the end of image response to ignore in ole test

psthBinTime = .001 ; % (sec) size of psthNi Bin step
psthSmoothTime = .02 ; % (sec) size of psthNi Bin (smooth sliding window)

ConcatPathNum = 1 ; %TEMP -SHOULD BE INPUT
ImPathNum = 1 ;
DsPathNum = 2 ;

% parameters
Color_list = {'k','r','b','g','c','y','k','r','b','g','c','y','k','r','b','g','c','y','k','r','b','g','c','y'} ; % order of colors for each 

saveDsIdsPath = ['/Users/jcafaro/Documents/AnalysisFigures/NatStimDs/DsIdsDb',num2str(DB),'ConcatPathNum',num2str(ConcatPathNum),'DsPathNum',num2str(DsPathNum)] ;
    
% load stimulus 
slashi = strfind(DataBlock(DB).ImSlide{ImPathNum},'/') ; % find the /
StimPath = [DataBlock(DB).ImSlide{ImPathNum}(1:slashi(end-1)),'stimuli/s',DataBlock(DB).ImSlide{ImPathNum}(end-1:end)] ;
load(StimPath) ;

% correct directions > 360
stimulus.directions(stimulus.directions>=360)=stimulus.directions(stimulus.directions>=360)-360 ;
stimulus.directions = sort(stimulus.directions) ;
stimulus.directionsShown(stimulus.directionsShown>=360)=stimulus.directionsShown(stimulus.directionsShown>=360)-360 ;

% load data
dataRun = load_data(DataBlock(DB).ImSlideConcat{ConcatPathNum}) ;
dataRun = load_neurons(dataRun) ;
dataRun.piece.array_id = 1551 ; % CAUTION
dataRun = load_ei(dataRun, 'all') ;
NumCells = length(dataRun.spikes) ; % number of cells

% ImSlide stimulus time
dataRunIm = load_data(DataBlock(DB).ImSlide{ImPathNum}) ;
dataRunIm = load_neurons(dataRunIm) ;
ImStimTime = dataRunIm.duration ;
DgStimTime = dataRun.duration-ImStimTime ;
clear dataRunIm ;

% get ei center of mass 
for cells = 1:length(dataRun.spikes) ; % for each cell
    EiCtr(cells,:) = get_ei_com(dataRun, dataRun.cell_ids(cells), numElectrodeLayers) ;
end

% get triggers from 
triggs = dataRun.triggers(dataRun.triggers>DgStimTime) ; % image slide

% stim params
if iscell(stimulus.image) ;
    NumImages = length(stimulus.image) ; % number of images shown
else
    NumImages = 1 ;
end
NumDirNi = size(stimulus.directionsShown,2) ; % number of directions
NumTrialsAttempted = size(stimulus.directionsShown,1) ;
NumStimShort = NumImages * NumDirNi * NumTrialsAttempted -length(triggs) ; % number of stimuli short

% look for single missing triggers and estimate
TriggOutlierFactor = 1.75 ;
if NumStimShort>0 ; 
    triggsFixed = triggs ;
    triggDiffTime = median(diff(triggs)) ;
    for t=2:length(triggs) ;
        if (triggs(t)-triggs(t-1)) > triggDiffTime*TriggOutlierFactor ; % if the trigger is much later than other triggers
            NewTrigg = (triggs(t)-triggs(t-1))/2 
            triggsFixed = [triggsFixed(1:t-1); NewTrigg; triggs(t:end)] ; % add new trigger
        end
    end
    triggs = triggsFixed ;
end

% take out trial sets that were not complete
NumCompleteTrials = floor(length(triggs)/(NumDirNi*NumImages)) ; % number of trials with all directions and all images
directionsShown = stimulus.directionsShown(1:NumCompleteTrials,:) ; % take out the entire unfinished trials
triggs = triggs(1:NumCompleteTrials*NumImages*NumDirNi) ; % triggers for stim included in completed trials
 
% expand directionsShown for mutliple images
directionsShownFull = [] ;
for t = 1:NumCompleteTrials ; % for each trial
    directionsShownFull = [directionsShownFull;repmat(directionsShown(t,:),[NumImages 1])];
end

% stim params
StimDuration = mean(diff(triggs)) ; % average stimulus duration
StimFrameRate = StimDuration/stimulus.num_frames ; % average frame rate sec/frame
dirShownFullTranspose = directionsShownFull' ;

% id ds cells in DsPath
Params.TimeBounds = [0,DgStimTime] ;
Params.DsPathNum = DsPathNum ;
try load(saveDsIdsPath) ; % if they are ds ids already saved
catch % if not find them

    DataBlock(DB).DsConcatPath = DataBlock(DB).ImSlideConcat{ConcatPathNum} ;
    ForIgor = DsCellFinder(DataBlock, DB, Params) ;

    save(saveDsIdsPath, 'ForIgor')
end

% On-Off cell from DsPath
cell_id = ForIgor.ds_id{2} ; 
DsTypeName = ForIgor.dsName{2} ; 
dsi = ForIgor.dsi{2} ;

for DsType=1:length(DsTypeName) ;
    DsCelli{DsType} = get_cell_indices(dataRun,cell_id{DsType}) ;
    DsCellTypeNum(DsType) = length(DsCelli{DsType}) ; 
end

% psthNi
psthTimeNi = [0:psthBinTime:StimDuration] ;
psthSmoothPnts = floor(psthSmoothTime/psthBinTime) ; % number of psthBins per smooth window
psthNi = nan(NumCells,NumDirNi,NumImages,NumCompleteTrials,length(psthTimeNi)) ;
for cells=1:NumCells ; % for each cell
    for st = 1:NumDirNi ; % for each stimulus direction
        sti = find(dirShownFullTranspose(:)==stimulus.directions(st)) ; % index of triggers for that direction
        for im = 1:NumImages ; % for each image
            ti = sti([im:NumImages:length(sti)]) ; % trigger index for that image and direction 
            
            for tr=1:(NumCompleteTrials) ; % for each trial
                spk = dataRun.spikes{cells}-triggs(ti(tr)) ;
                spk = spk(spk>=0 & spk<=StimDuration) ;

                spikeTrain = zeros(1,length(psthTimeNi)) ; % empty spike train
                spkPnts = 1+floor(spk/psthBinTime) ; % spike points 
                spikeTrain(spkPnts) = 1 ; % spike train
                psthNi(cells,st,im,tr,:) = smooth(spikeTrain,psthSmoothPnts) ; % psthNi (cells,dir,image,trial,time)
            end
        end
    end
end
psthNi = psthNi/psthBinTime ; % change to hz

%% mean and variance across trials or time/image 

% points and permutations
Tpnts = [StimOnOffset/psthBinTime:(length(psthTimeNi)-StimOffTime/psthBinTime)] ; % selected time points
psthNiSel = psthNi(:,:,:,:,Tpnts) ; % select time subsets - psthNiSel(cells,dir,image,trial,time)
psthNiPerm = permute(psthNiSel,[1,2,4,3,5]) ; % psthNiSel(cells,dir,image,trial,time)-->psthNiPerm(cells,dir,trial,image,time)
psthNiPermRe = reshape(psthNiPerm,[NumCells,NumDirNi,NumCompleteTrials,(NumImages*length(Tpnts))]) ; % psthNiPermRe(cells,dir,trial,image/time) ;

% across trials
psthNi_varAcrossTrials = squeeze(var(psthNiPermRe,[],3)) ; % (cells,dir,image/time)
psthNi_meanAcrossTrials = squeeze(mean(psthNiPermRe,3)) ; % (cells,dir,images/time)

psthNi_varAcrossTrials_mean = mean(psthNi_varAcrossTrials,3) ; % (cells,dir)
psthNi_meanAcrossTrials_mean = mean(psthNi_meanAcrossTrials,3) ; % (cells,dir)

psthNi_meanAcrossTrials_var = var(psthNi_meanAcrossTrials,[],3) ; % (cells,dir)

% across time/image
psthNiPermRe_varAcrossTimeAndImages = squeeze(var(psthNiPermRe,[],4)) ; % (cells,dir,trial)
psthNiPermRe_meanAcrossTimeAndImages = squeeze(mean(psthNiPermRe,4)) ; % (cells,dir,trial)

psthNiPermRe_varAcrossTimeAndImages_mean = mean(psthNiPermRe_varAcrossTimeAndImages,3) ; % (cells,dir)
psthNiPermRe_meanAcrossTimeAndImages_mean = mean(psthNiPermRe_meanAcrossTimeAndImages,3) ; % (cells,dir)

psthNiPermRe_meanAcrossTimeAndImages_var = var(psthNiPermRe_meanAcrossTimeAndImages,[],3) ; % (cells,dir)

%% structure of noise

% pca
[evecs,evals] = eig(cvMat) ;

% pdfs for example cell
egi = DsCelli{2}(5) ; % example cell i 
ExamplePdfX = [0:6:60] ; 
ExampleTc = psthNi_meanAcrossTrials_mean(egi,:) ; % eg tuning curve (mean responses)
for st = 1:NumDirNi ; % for each stimulus direction
    [ExamplePdf(st,:)] = hist(squeeze(psthNi_meanAcrossTrials(egi,st,:)),ExamplePdfX) ; % (cells,dir,images/time)
    ExamplePdf(st,:) = ExamplePdf(st,:)/sum(ExamplePdf(st,:)) ; % make pdf
end
    
%% fisher information - pairs
FlatFlag = 1 ;

% points and permutations (if not run above)
Tpnts = [StimOnOffset/psthBinTime:(length(psthTimeNi)-StimOffTime/psthBinTime)] ; % selected time points
psthNiSel = psthNi(:,:,:,:,Tpnts) ; % select time subsets - psthNiSel(cells,dir,image,trial,time)
psthNiPerm = permute(psthNiSel,[1,2,4,3,5]) ; % psthNiSel(cells,dir,image,trial,time)-->psthNiPerm(cells,dir,trial,image,time)
psthNiPermRe = reshape(psthNiPerm,[NumCells,NumDirNi,NumCompleteTrials,(NumImages*length(Tpnts))]) ; % psthNiPermRe(cells,dir,trial,image/time) ;

psthNi_meanAcrossTrials = squeeze(mean(psthNiPermRe,3)) ; % (cells,dir,images/time)
psthNi_meanAcrossTrials_mean = mean(psthNi_meanAcrossTrials,3) ; % (cells,dir)
psthNi_meanAcrossTrials_mean_circ = [psthNi_meanAcrossTrials_mean,psthNi_meanAcrossTrials_mean(:,1)] ; % circular

% info
Fi = zeros(NumCells) ; % prep mat
FiInd = zeros(NumCells) ;
cvMat = zeros(NumCells) ;
DistMat = nan(NumCells) ;
PairMat = nan(NumCells) ;

for c1=1:NumCells-1 ; % for each cell
    for c2=(c1+1):NumCells ; % for each cell
        Ipair(N) = get_MI_from_pair_resp(Rpair,[1:ns],1,[1:nt]) ;
        
        DistMat(c1,c2) = sqrt(sum((EiCtr(c1,:)-EiCtr(c2,:)).^2)) ; % (um) distance between cells
        PairMat(c1,c2) = ismember(c1,cell2mat(DsCelli)) + ismember(c2,cell2mat(DsCelli)) ; % 0=non-non, 1=ds-non, 2=ds-ds

        for st=1:NumDirNi ; % for each direction
            diff1 = psthNi_meanAcrossTrials_mean_circ(c1,st)-psthNi_meanAcrossTrials_mean_circ(c1,st+1) ;
            diff2 = psthNi_meanAcrossTrials_mean_circ(c2,st)-psthNi_meanAcrossTrials_mean_circ(c2,st+1) ;
             
            if FlatFlag ;
                diff1 = diff1*ismember(c1,cell2mat(DsCelli)) ; % artificially set change to zero
                diff2 = diff2*ismember(c2,cell2mat(DsCelli)) ;
            end
                 
            cv = cov(squeeze(psthNi_meanAcrossTrials(c1,st,:)),squeeze(psthNi_meanAcrossTrials(c2,st,:))) ;
            cvInd = cv.*eye(2) ; % get rid of cov keep variance
            
            if c2==c1+1 ; % prevents repition of variance measure
                cvMat(c2,c2) = cvMat(c2,c2) + cv(2,2) ; % variance
                if c1==1 ; % prevents repition of variance measure
                    cvMat(c1,c1) = cvMat(c1,c1) + cv(1,1) ; % variance
                end
            end

            cvMat(c1,c2) = cvMat(c1,c2) + cv(1,2) ; % covariance
            
            Fi(c1,c2) = Fi(c1,c2) + [diff1,diff2]*inv(cv)*[diff1,diff2]' ; % fisher info
            FiInd(c1,c2) = FiInd(c1,c2) + [diff1,diff2]*inv(cvInd)*[diff1,diff2]' ; % fisher info if the cells were independent
        end
    end
end
cvMat = (cvMat+triu(cvMat,1)')/NumDirNi ; % average cov mat within direction (ie distractor noise cov mat)     
Ishuff = Fi - FiInd ; % I shuff
IshuffFrac = Ishuff./(Fi + FiInd) ; % info change relative to total info

corrMat = corrcov(cvMat) ; % corr coefs

%% Fisher info - population

% points and permutations (if not run above)
Tpnts = [StimOnOffset/psthBinTime:(length(psthTimeNi)-StimOffTime/psthBinTime)] ; % selected time points
psthNiSel = psthNi(:,:,:,:,Tpnts) ; % select time subsets - psthNiSel(cells,dir,image,trial,time)
psthNiPerm = permute(psthNiSel,[1,2,4,3,5]) ; % psthNiSel(cells,dir,image,trial,time)-->psthNiPerm(cells,dir,trial,image,time)
psthNiPermRe = reshape(psthNiPerm,[NumCells,NumDirNi,NumCompleteTrials,(NumImages*length(Tpnts))]) ; % psthNiPermRe(cells,dir,trial,image/time) ;

psthNi_meanAcrossTrials = squeeze(mean(psthNiPermRe,3)) ; % (cells,dir,images/time)
psthNi_meanAcrossTrials_mean = mean(psthNi_meanAcrossTrials,3) ; % (cells,dir)
psthNi_meanAcrossTrials_mean_circ = [psthNi_meanAcrossTrials_mean,psthNi_meanAcrossTrials_mean(:,1)] ; % circular

% d^2 (Fisher Info - Averbeck 2006)
NonDsi = setdiff([1:NumCells],cell2mat(DsCelli)) ;
cvMat = zeros(NumCells) ;
for st=1:NumDirNi ; % for each direction
    diffs = psthNi_meanAcrossTrials_mean_circ(:,st)-psthNi_meanAcrossTrials_mean_circ(:,st+1) ;
    diffsDsOnly = psthNi_meanAcrossTrials_mean_circ(cell2mat(DsCelli),st)-psthNi_meanAcrossTrials_mean_circ(cell2mat(DsCelli),st+1) ;
    diffsNonDsOnly = psthNi_meanAcrossTrials_mean_circ(NonDsi,st)-psthNi_meanAcrossTrials_mean_circ(NonDsi,st+1) ;
    
    CovData = squeeze(psthNi_meanAcrossTrials(:,st,:)) ; % data for cov mat 
    CovDataDsOnly = squeeze(psthNi_meanAcrossTrials(cell2mat(DsCelli),st,:)) ; % data for cov mat 
    CovDataNonDsOnly = squeeze(psthNi_meanAcrossTrials(NonDsi,st,:)) ; % data for cov mat 
    
    cv = cov(CovData') ;
    cvInd = cv.*eye(NumCells) ; % get rid of cov keep variance
    cvMat = cvMat+cv/NumDirNi ; % for posterity (should be the same as above)
    cvDsOnly = cov(CovDataDsOnly') ; % DS cells only
    cvNonDsOnly = cov(CovDataNonDsOnly') ; % DS cells only

    FiPop(st) = diffs'*inv(cv)*diffs ; % fisher info
    FiPopInd(st) = diffs'*inv(cvInd)*diffs ; % fisher info if the cells were independent 
    FiPopDiag(st) = (FiPopInd(st)^2)/(diffs'*inv(cvInd)*cv*inv(cvInd)*diffs) ; % fisher info if trained on shuffle
    FiPopDsOnly(st) = diffsDsOnly'*inv(cvDsOnly)*diffsDsOnly ; % fisher info
    FiPopNonDsOnly(st) = diffsNonDsOnly'*inv(cvNonDsOnly)*diffsNonDsOnly ; % fisher info
end
        
PcPop = 1-erfc(sqrt(FiPop)/2) ;    
PcPopInd = 1-erfc(sqrt(FiPopInd)/2) ;
PcPopDiag = 1-erfc(sqrt(FiPopDiag)/2) ;
PcPopDsOnly = 1-erfc(sqrt(FiPopDsOnly)/2) ;
PcPopNonDsOnly = 1-erfc(sqrt(FiPopNonDsOnly)/2) ;

%% Idiag - OLE for all stim directions

psthNi_meanAcrossTrials_Perm = permute(psthNi_meanAcrossTrials,[1,3,2]) ; % (cells,image/time,dir)

TrainPnts = [1:size(psthNi_meanAcrossTrials_Perm,2)/2] ;
TestPnts = [size(psthNi_meanAcrossTrials_Perm,2)/2+1:size(psthNi_meanAcrossTrials_Perm,2)] ;

TrainData = reshape(psthNi_meanAcrossTrials_Perm(:,TrainPnts,:) ,[NumCells,NumDirNi*length(TrainPnts)]) ; % (cells,dir/images/time) 
TestData = reshape(psthNi_meanAcrossTrials_Perm(:,TestPnts,:) ,[NumCells,NumDirNi*length(TestPnts)]) ;
for cells=1:NumCells ;
    TrainDataShuff(cells,:) = reshape(psthNi_meanAcrossTrials_Perm(cells,TrainPnts(randperm(length(TrainPnts))),:) ,[1,NumDirNi*length(TrainPnts)]) ; %
    TestDataShuff(cells,:) = reshape(psthNi_meanAcrossTrials_Perm(cells,TestPnts(randperm(length(TestPnts))),:) ,[1,NumDirNi*length(TestPnts)]) ; % 
end
    
DirData(1,:,:) = repmat(stimulus.directions,[size(psthNi_meanAcrossTrials_Perm,2),1]) ;
TrainDirVector = reshape(DirData(1,TrainPnts,:) ,[1,NumDirNi*length(TrainPnts)]) ; 
TestDirVector = reshape(DirData(1,TestPnts,:) ,[1,NumDirNi*length(TestPnts)]) ; 

[OleWeights, OlePolarWeights] = OleFinder(TrainDirVector,TrainData') ;
[DirEstimate,DirEstimate_mag]= OleTester(OleWeights,TestData') ;

for tp=1:length(DirEstimate) ;% error
    Error(tp) = acuteAngle(DirEstimate(tp),TestDirVector(tp)) ;
end
 
% train shuffle, test actual
[OleWeightsShuff, OlePolarWeightsShuff] = OleFinder(TrainDirVector,TrainDataShuff') ;

[DirEstimateShuff,DirEstimateShuff_mag]= OleTester(OleWeightsShuff,TestData') ;

for tp=1:length(DirEstimate) ;% error
    ErrorShuff(tp) = acuteAngle(DirEstimateShuff(tp),TestDirVector(tp)) ;
end

% train shuffle, test shuffle
[DirEstimateShuff2,DirEstimateShuff2_mag]= OleTester(OleWeightsShuff,TestDataShuff') ;

for tp=1:length(DirEstimate) ;% error
    ErrorShuff2(tp) = acuteAngle(DirEstimateShuff2(tp),TestDirVector(tp)) ;
end

% Error as function of estimate mag
cnt = 1 ;
for DmThresh = [0:.1:5] ;
    ErrorVsMag(cnt) = median(Error(DirEstimate_mag>DmThresh)) ;
    ErrorVsMagShuff(cnt) = median(ErrorShuff(DirEstimateShuff_mag>DmThresh)) ;
    ErrorVsMagShuff2(cnt) = median(ErrorShuff2(DirEstimateShuff2_mag>DmThresh)) ;
    cnt = cnt+1 ;
end
    
%% modified OLE

alpha = 0.0005 ;
beta = 0.000005 ;

TempM = (cvMat*TestData)' ;
ModAlpha = -alpha*TempM*OleWeights(1:end-1,:) ;

TempM = (cvMat*(TestData))' ;
pX = TestData.*repmat(OleWeights(1:end-1,1),[1,size(TestData,2)]) ;
for pnt=1:size(pX,2) ; % for each point
    ModBetaX(pnt) = beta*TempM(pnt,:)*pX(:,pnt) ;
end

pY = TestData.*repmat(OleWeights(1:end-1,2),[1,size(TestData,2)]) ;
for pnt=1:size(pX,2) ; % for each point
    ModBetaY(pnt) = beta*TempM(pnt,:)*pY(:,pnt) ;
end

ModBeta = [ModBetaX',ModBetaY'] ;

CartMod = ModAlpha + ModBeta ;

% polar estimate
PolarMod = atan2d(CartMod(:,2), CartMod(:,1)) ; % four quadrant inverse tangent
PolarMod(PolarMod<0) = PolarMod(PolarMod<0)+360 ; % no negatives
PolarModMag = sqrt(CartMod(:,2).^2 + CartMod(:,1).^2) ;

for imtb = 1:length(DirEstimate) ; % for each image time bin
    Temp = PolarVectorAddition([DirEstimate(imtb),DirEstimate_mag(imtb); PolarMod(imtb),PolarModMag(imtb)]) ; 
    DirEstimateMod(imtb) = Temp(1) ;
    DirEstimateMod_mag(imtb) = Temp(2) ;
end


%% assess info lost for each cell

FiPopExclude1 = nan(NumCells,NumDirNi) ; % prep mat
Fi1 = nan(NumCells,NumDirNi) ; % prep mat

for cells = 1:NumCells ; % for each cell
    TempCells = setdiff([1:NumCells],cells) ;
    for st=1:NumDirNi ; % for each direction
        diffsMod = psthNi_meanAcrossTrials_mean_circ(TempCells,st)-psthNi_meanAcrossTrials_mean_circ(TempCells,st+1) ;
        
        CovData = squeeze(psthNi_meanAcrossTrials(TempCells,st,:)) ; % data for cov mat 

        cv = cov(CovData') ;

        FiPopExclude1(cells,st) = diffsMod'*inv(cv)*diffsMod ; % 
        
        diffs1 = psthNi_meanAcrossTrials_mean_circ(cells,st)-psthNi_meanAcrossTrials_mean_circ(cells,st+1) ;
        Var1 = var(squeeze(psthNi_meanAcrossTrials(cells,st,:))) ;
        Fi1(cells,st) = diffs1*inv(Var1)*diffs1 ; %
    end
end




%% figures

figure % psth
cells = 85 ;
st = 1 ;
im = 1 ;
for st=1:NumDirNi ;
    subplot(NumDirNi,1,st)
    PlotOffset = 0 ;
    for tr=1:NumCompleteTrials;
        plot(psthTimeNi,squeeze(psthNi(cells,st,im,tr,:))+PlotOffset)
        hold on
        PlotOffset = PlotOffset+max(squeeze(psthNi(cells,st,im,tr,:)))+1 ;
    end
end

figure % mean vs variance 
for DsType = 1:length(DsCelli) ;
    for cells = 1:length(DsCelli{DsType}) ;
        subplot(1,3,1) % mean vs var - CAUTION could be too many points for ploter
        plot(psthNiPermRe_meanAcrossTimeAndImages(DsCelli{DsType}(cells),:),...
            psthNiPermRe_varAcrossTimeAndImages(DsCelli{DsType}(cells),:),'r.')
        hold on
        plot(psthNi_meanAcrossTrials(DsCelli{DsType}(cells),:),...
            psthNi_varAcrossTrials(DsCelli{DsType}(cells),:),'k.')
        xlabel('mean')
        ylabel('variance')
        
        subplot(1,3,2) % mean vs var - averaged within direction
        plot(psthNiPermRe_meanAcrossTimeAndImages_mean(DsCelli{DsType}(cells),:),...
            psthNiPermRe_varAcrossTimeAndImages_mean(DsCelli{DsType}(cells),:),'r.')
        hold on
        plot(psthNi_meanAcrossTrials_mean(DsCelli{DsType}(cells),:),...
            psthNi_varAcrossTrials_mean(DsCelli{DsType}(cells),:),'k.')
        xlabel('mean')
        ylabel('variance')
        
        subplot(1,3,3) % mean vs var - mean and variance of averages
        plot(psthNiPermRe_meanAcrossTimeAndImages_mean(DsCelli{DsType}(cells),:),...
            psthNiPermRe_meanAcrossTimeAndImages_var(DsCelli{DsType}(cells),:),'k.')
        hold on
        plot(psthNi_meanAcrossTrials_mean(DsCelli{DsType}(cells),:),...
            psthNi_meanAcrossTrials_var(DsCelli{DsType}(cells),:),'r.')
        xlabel('mean')
        ylabel('variance')
    end
end
plot([0,20],[0,20],'g')
xlabel('mean')
ylabel('variance')

figure % example cell and Pdf
subplot(2,1,1)
plot(ExampleTc)

subplot(2,1,2)
for st = 1:NumDirNi ; % for each stimulus direction
    plot(ExamplePdfX,ExamplePdf(st,:))
    hold on
end

figure % each trial or time/image point 
for DsType = 1:length(DsCelli) ;
    for cells = 1:length(DsCelli{DsType}) ;
        subplot(2,1,1)
        for trial = 1:NumCompleteTrials ;
            plot(squeeze(psthNi_avAcrossTimeAndImages(DsCelli{DsType}(cells),:,trial)),'*')
            hold on
        end
        hold off
        
        subplot(2,1,2)
        for pnt = 1:length(Tpnts) ; % for each time/image point 
            plot(squeeze(psthNi_avAcrossTrials_comb(DsCelli{DsType}(cells),:,pnt)),'.')
            hold on
        end
        hold off
        pause
    end   
end

figure % Ishuff vs dist -pairs
subplot(3,1,1) % ds-non pair
plot(DistMat(PairMat==1),Ishuff((PairMat==1)),'*')

subplot(3,1,2) % ds-ds pair
plot(DistMat(PairMat==2),Ishuff((PairMat==2)),'*')

subplot(3,1,3) % non-non pair
plot(DistMat(PairMat==0),Ishuff((PairMat==0)),'*')

figure % Ishuff vs corrcoef -pairs
subplot(3,1,1) % ds-non pair
plot(corrMat(PairMat==1),Ishuff(PairMat==1)./Fi((PairMat==1)),'*')
hold on
[X,Y,Ysem] = NonLinFilterFinder(corrMat(PairMat==1),Ishuff(PairMat==1)./Fi((PairMat==1)),0.05) ;
plot(X,Y)
xlabel('corr coef')
ylabel('Ishuff')

subplot(3,1,2) % ds-ds pair
plot(corrMat(PairMat==2),Ishuff((PairMat==2))./Fi((PairMat==2)),'*')
hold on
[X,Y,Ysem] = NonLinFilterFinder(corrMat(PairMat==2),Ishuff(PairMat==2)./Fi((PairMat==2)),0.05) ;
plot(X,Y)
xlabel('corr coef')
ylabel('Ishuff')

subplot(3,1,3) % non-non pair
plot(corrMat(PairMat==0),Ishuff((PairMat==0))./Fi((PairMat==0)),'*')
hold on
[X,Y,Ysem] = NonLinFilterFinder(corrMat(PairMat==0),Ishuff(PairMat==0)./Fi((PairMat==0)),0.05) ;
plot(X,Y)
xlabel('corr coef')
ylabel('Ishuff')

figure % dimensionality of distractor noise
subplot(2,1,1)
plot(diag(evals)/sum(diag(evals)),'*')
xlabel('pc')
ylabel('fraction variance')

subplot(2,1,2)
plot(evecs(:,end))
hold on
plot(diag(cvMat)/max(diag(cvMat)))
xlabel('cell')
ylabel('weight or variance')

figure % Percent correct for binary classifier
plot(PcPop,'b')
hold on
plot(PcPopInd,'c')
plot(PcPopDiag,'y')
plot(PcPopDsOnly,'r')
plot(PcPopNonDsOnly,'r:')
xlabel('direction set')
ylabel('percent correct')
legend('full- encoded','shuffled-encoded','trained-shuff full-decoded','Ds only','NonDsOnly')

figure % Idiag - OLE
subplot(4,1,1)
plot(DirEstimateShuff,'b')
hold on
plot(DirEstimate,'r')
plot(DirEstimateShuff2,'y')
plot(TestDirVector,'k')
xlabel('observation')
ylabel('direction (deg)')

subplot(4,1,2) % error hist
plot([0:5:180],cumsum(hist(Error,[0:5:180]))/sum(hist(Error,[0:5:180])),'r')
hold on
plot([0:5:180],cumsum(hist(ErrorShuff,[0:5:180]))/sum(hist(ErrorShuff,[0:5:180])),'b')
plot([0:5:180],cumsum(hist(ErrorShuff2,[0:5:180]))/sum(hist(ErrorShuff2,[0:5:180])),'y')
xlabel('error')
ylabel('fraction observed')

subplot(4,1,3) % estimate mag 
plot([0:.1:5],cumsum(hist(DirEstimate_mag,[0:.1:5]))/sum(hist(DirEstimate_mag,[0:.1:5])),'r')
hold on
plot([0:.1:5],cumsum(hist(DirEstimateShuff_mag,[0:.1:5]))/sum(hist(DirEstimateShuff_mag,[0:.1:5])),'b')
plot([0:.1:5],cumsum(hist(DirEstimateShuff2_mag,[0:.1:5]))/sum(hist(DirEstimateShuff2_mag,[0:.1:5])),'y')
xlabel('estimate magnitude')
ylabel('fraction observed')

subplot(4,1,4) % median error as function of estimate mag
plot([0:.1:5],ErrorVsMag,'r')
hold on
plot([0:.1:5],ErrorVsMagShuff,'b')
plot([0:.1:5],ErrorVsMagShuff2,'y')

% Fisher with cells excluded
figure
subplot(2,1,1) % DS and nonDS
plot((sum(FiPop)-sum(FiPopExclude1,2))/sum(FiPop),'*')
hold on
plot(cell2mat(DsCelli),(sum(FiPop)-sum(FiPopExclude1(cell2mat(DsCelli),:),2))/sum(FiPop),'ro')
xlabel('cell')
ylabel('Fraction FI lost')

subplot(2,1,2) % 
plot(sum(Fi1,2),(sum(FiPop)-sum(FiPopExclude1,2))/sum(FiPop),'*')
hold on
plot(sum(Fi1(cell2mat(DsCelli),:),2),(sum(FiPop)-sum(FiPopExclude1(cell2mat(DsCelli),:),2))/sum(FiPop),'ro')
xlabel('FI solo')
ylabel('Fraction Fisher lost')




function ForIgor = RepStimClusterAnalysis(DataBlock, DB, Params) 

% JC 6/15/2017
%This function will compare clustering of stimulus (may be repeated) and response

mapEiFlag = false ;
SquareGridSizeUpperBound = 300 ;% (um) ~size of boxes on one side of a square grid dividing up array space
numElectrodeLayers = 2 ; % number of electrode layers used to calc ei center
RepStimNum = 1 ; % select repstim set
BwPathNum = 1 ;
numTrigPerRep = 3 ; % (integer) number of triggs per stimulus repeat
PsthBinTime = 0.01 ; % (sec) 
BinSearchNumber = 2 ; % (integer) number of psth bins that can be compared 
StimFrameNumber = 300 ; % (frames) number of frames in repeated stimulus
StimFrameStimNumber = 10 ; % (frames) number of frames that RGCs can integrate over
StimPath = '/Volumes/lab/Documents/Movies/CatCam/cat_mean117_sd62_0to255' ; 
numClusters=[1,5,10,25,50,100] ; % k in k-means

% load data
dataRun = load_data(DataBlock(DB).RepStim{RepStimNum}) ;
dataRun = load_neurons(dataRun) ;
dataRun = load_ei(dataRun,'all') ;

% load BW data with identified neurons
dataRunBw = load_data(DataBlock(DB).BwPath{BwPathNum}) ;
dataRunBw = load_neurons(dataRunBw) ;
dataRunBw = load_params(dataRunBw,'cell_type_depth', 5) ;
     
% identified cell types

for a = 1:length(dataRunBw.cell_types) ;
    celltypes{a} = dataRunBw.cell_types{a}.name ;
end

UniqueCellTypes = unique(celltypes) ;
if isempty(UniqueCellTypes{1}) ;
    UniqueCellTypes = UniqueCellTypes(2:end) ;
end

% map dim flash onto bw data
if mapEiFlag ; % if using map ei cells

    % map using electrical images
    cell_list_map = map_ei(dataRunBw, dataRun) ;

    % cells ids in slave for each UniqueCellType set in master data
    for uc = 1:length(UniqueCellTypes) ;
        Masteri{uc} = get_cell_indices(dataRunBw, UniqueCellTypes{uc}) ;
        cell_ids{uc} = cell2mat(cell_list_map(Masteri{uc})) ;
    end
else % if not using map ei
    for uc = 1:length(UniqueCellTypes) ;
        cell_ids{uc} = intersect(dataRun.cell_ids, get_cell_ids(dataRunBw,UniqueCellTypes{uc})) ;
        cell_i{uc} = get_cell_indices(dataRun,cell_ids{uc}) ;
    end
end

%% find joint histogram from clusters - from entire array and stim area

% response
TrialStartTimes = dataRun.triggers([1:numTrigPerRep:length(dataRun.triggers)]) ;
TrialDuration = min(diff(TrialStartTimes)) ;
trialNumber = length(TrialStartTimes) ; % number of trials
numCells = length(dataRun.spikes) ; % number of cells
PsthTime = [0:PsthBinTime:TrialDuration] ; % (s)

r1 = cell(1,trialNumber) ; % prep response vector
r2 = cell(1,trialNumber) ; % prep response vector

for tr = 1:trialNumber ; % for each trial
    Bins = PsthTime+TrialStartTimes(tr) ; % time bins for trial 
    for c = 1:numCells ; % for each cell
        r{tr}(:,c)=histcounts(dataRun.spikes{c},Bins) ; % number of spikes each bin
    end
end

% stimulus
load(StimPath) ; % movie
Stim = mov(:,:,1:StimFrameNumber) ; % section of movie shown
[movPCs,movPrj,PcVal] = MovPCA(Stim,StimFrameStimNumber) ; % PCA on movie stimulus

movPrjWtd = movPrj.*repmat(PcVal',size(movPrj,1),1) ; % wieght Prj values by fraction variance 

StimTime = [StimFrameStimNumber:length(movPrj)+StimFrameStimNumber-1]...
    *(PsthTime(end)/(length(movPrj)+StimFrameStimNumber-1)) ; % time of stimulus FIXXXX!!!!!!!!!!

movPrjWtd_interp = interp1(StimTime,movPrjWtd,PsthTime) ; % interpolate weighted stimulus projection values

% k-means clustering 
rnd=0 ;
for k = numClusters ; % for each k in k-means
    rnd=rnd+1 ;
    StimClusters(:,rnd) = kmeans(movPrjWtd_interp,k) ;
    ResponseClusters(:,rnd) = kmeans(cell2mat(r'),k) ;
end 

% joint distribution (column-->count of response cluster given a stim cluster)
RbinIds = repmat([1:length(PsthTime)-1],1,trialNumber) ; % time bins of ResponseClusters (because r was concatinated over trials)

for k=1:length(numClusters) ; % for each k in k-means
    JointDist{k}= zeros(numClusters(k)) ; % initialize Joint Distribution
    for clstr = 1:numClusters(k) ; % for each cluster
        clstri = find(StimClusters(:,k) == clstr) ; % time bins within stim cluster 
        for pnt = 1:length(clstri) ; % for each time bin in the stim cluster
            JointDist{k}(:,clstr) = JointDist{k}(:,clstr) + ...
                hist(ResponseClusters(RbinIds==clstri(pnt),k),[1:numClusters(k)])' ; % add the histogram 
        end
    end
end
            
% precent correct for optimal ML discrimitor
for k=1:length(numClusters) ; % for each k in k-means
    PercentCorrect(k) = sum(max(JointDist{k},[],2))/sum(JointDist{k}(:)) ;
end


%% find joint histogram from clusters - in small groups

% calculate ei centers
for c = 1:length(dataRun.spikes) ; % for each cell
    eiCnt(c,:) = get_ei_com(dataRun, dataRun.cell_ids(c), numElectrodeLayers) ;
end

% find groups of cells in each grid bin

SquareGridNum_X = ceil(range(dataRun.ei.array_bounds_x)/SquareGridSizeUpperBound) ; % bin number
SquareGridNum_Y = ceil(range(dataRun.ei.array_bounds_y)/SquareGridSizeUpperBound) ;

SquareGridSize_X = range(dataRun.ei.array_bounds_x)/SquareGridNum_X ; % bin size
SquareGridSize_Y = range(dataRun.ei.array_bounds_y)/SquareGridNum_Y ;

r=1 ;
for x=1:SquareGridNum_X ;
    for y = 1:SquareGridNum_Y ;
        xbounds = [(x-1)*SquareGridSize_X+dataRun.ei.array_bounds_x(1),...
            x*SquareGridSize_X+dataRun.ei.array_bounds_x(1)] 
        ybounds = [(y-1)*SquareGridSize_Y+dataRun.ei.array_bounds_y(1),...
            y*SquareGridSize_Y+dataRun.ei.array_bounds_y(1)] 
        
        xi = find(eiCnt(:,1)>=xbounds(1) & eiCnt(:,1)<xbounds(2)) ;
        yi = find(eiCnt(:,2)>=ybounds(1) & eiCnt(:,2)<ybounds(2)) ;
        Group_i{r} = intersect(xi,yi) ; % indicies of cells within the box
        r=r+1 ;
    end
end
   

for g=1:length(Group_i) ; % for each group
    PopDist_groups_struct = PopDistFinder(dataRun.spikes(Group_i{g}), dataRun.triggers([1:numTrigPerRep:length(dataRun.triggers)]),...
        'PsthBinTime', PsthBinTime,'BinSearchNumber', BinSearchNumber, 'TrialNumberMax',3) ;
    PopDist_groups{g} = mean(PopDist_groups_struct.AcrossStimMinusAcrossTrials,1) ;    

    for uc = 1:length(UniqueCellTypes) ; % for each unique cell type
        if ~isempty(cell_i{uc}) ; % if there are cells
            gi = intersect(Group_i{g},cell_i{uc}) ; % type of cell within small group
            if ~isempty(gi) ;
                PopDist_groups_withOnly_struct = PopDistFinder(dataRun.spikes(gi), dataRun.triggers([1:numTrigPerRep:length(dataRun.triggers)]),...
                    'PsthBinTime', PsthBinTime,'BinSearchNumber', BinSearchNumber, 'TrialNumberMax',3) ;

                ci = setdiff(Group_i{g},gi) ; % indicies of all other cells
                PopDist_groups_without_struct = PopDistFinder(dataRun.spikes(ci), dataRun.triggers([1:numTrigPerRep:length(dataRun.triggers)]),...
                    'PsthBinTime', PsthBinTime,'BinSearchNumber', BinSearchNumber, 'TrialNumberMax',3) ;

                PopDist_groups_withOnly{g}{uc} = mean(PopDist_groups_withOnly_struct.AcrossStimMinusAcrossTrials,1) ;
                PopDist_groups_without{g}{uc} = mean(PopDist_groups_without_struct.AcrossStimMinusAcrossTrials,1) ;

                Necessity_groups{g}(uc,:) = (PopDist_groups{g}-PopDist_groups_without{g}{uc})./PopDist_groups{g} ; 
                Sufficiency_groups{g}(uc,:) = PopDist_groups_withOnly{g}{uc}./PopDist_groups{g};
            else
                PopDist_groups_withOnly{g}{uc} = nan ;
                PopDist_groups_without{g}{uc} = nan ;

                Necessity_groups{g}(uc,:) = nan(1,length(Necessity_all{1})) ; 
                Sufficiency_groups{g}(uc,:) = nan(1,length(Necessity_all{1})) ; 
            end
                
        end
    end
end
    

% discriminability of stimulus

load(StimPath) ; % movie
Stim = mov(:,:,1:StimFrameNumber) ; % section of movie shown
[movPCs,movPrj,PcVal] = MovPCA(Stim,StimFrameStimNumber) ; % PCA on movie stimulus

movPrjWtd = movPrj.*repmat(PcVal',size(movPrj,1),1) ; % wieght Prj values by fraction variance 

for f1=1:size(movPrj,1) ;
    StimDist_all(f1) = 0 ; % initialize
    for f2=1:size(movPrj,1) ;
        distTemp = sqrt(sum((movPrjWtd(f1,:)-movPrjWtd(f2,:)).^2)) ; % eucleidean distance between two stim 
        StimDist_all(f1) = StimDist_all(f1)+distTemp ; % sum of distances between movPrj(f1) and all other movPrj
    end
end

StimTime = [StimFrameStimNumber:length(StimDist_all)+StimFrameStimNumber-1]...
    *(PsthTime(end)/(length(StimDist_all)+StimFrameStimNumber-1)) ; % time of stimulus

movPrjWtd_interp = interp1(StimTime,movPrjWtd,PsthTime) ; % interpolate weighted stimulus projection values

% k-means clustering 

for k=5:100 ; % for each k in k-means
    StimClusters(:,k) = kmeans(movPrjWtd_interp,k) ;
    ResponseClusters(:,k) = kmeans(movPrjWtd_interp,k) ;
end
    
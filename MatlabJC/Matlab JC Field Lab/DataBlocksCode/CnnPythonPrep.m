function Temp = CnnPythonPrep(DataBlock, DB, Params) 

% load BW data block
dataRun = load_data(DataBlock(DB).BwPath{Params.BwPathNum}) ;
dataRun = load_neurons(dataRun) ;
%dataRun = load_ei(dataRun, 'all') ;
dataRun = load_params(dataRun,'cell_type_depth', 5) ;
dataRun = load_sta(dataRun) ; % only necessary to get trf_time
dataRun = get_sta_summaries(dataRun,'all') ;    

dataRun.stimulus.monitor_refresh = 60.35 ; % the monitor_refresh rate is not accurate

% get movie
dataRun = load_java_movie(dataRun,DataBlock(DB).BwMoviePath{Params.BwPathNum}) ;
movTemp = get_movie(DataBlock(DB).BWMoviePath{Params.BwPathNum}, dataRun.triggers, size(dataRun.stimulus.java_movie)) ;
mov(:,:,:) = movTemp(:,:,1,1:end-100) ; % get rid of rgb for achromatic stim and clip set frames after last trigger 

% make time bins consistent with estimated frame bins
FramesPerTrig = 100 ; % frames per triger
trigdiff = diff(dataRun.triggers) ;
PsthTime = dataRun.triggers(1) ;
for t=1:length(trigdiff) ; % for each trigger interval
    PsthTime = [PsthTime,PsthTime(end)+[1:FramesPerTrig]*trigdiff(t)/FramesPerTrig] ;
end
  
% make psth
NumCells = length(dataRun.spikes) ;
Psth = nans(NumCells,length(PsthTime)-1) ;
for c = 1:NumCells ; % for each cell
    Psth(c,:) = histcounts(dataRun.spikes{c},PsthTime)./diff(PsthTime) ; % spikes per sec
end
Psth = [zeros(NumCells,1),Psth] ; 

% arrange Psth by cell type
for a = 1:length(dataRun.cell_types) ;
    celltypes{a} = dataRun.cell_types{a}.name ;
end

UniqueCellTypes = unique(celltypes) ;
if isempty(UniqueCellTypes{1}) ;
    UniqueCellTypes = UniqueCellTypes(2:end) ;
end

Psth_bytypes = nans(size(Psth)) ; % prep mat
c=1 ;
for uc = 1:length(UniqueCellTypes) ; 
    Celltypei{uc} = get_cell_indices(dataRun, UniqueCellTypes{uc}) ;
    Psth_bytypei{uc} = c:c+length(Celltypei{uc})-1 ; % indecies in Psth_bytype mat
    Psth_bytypes(Psth_bytypei{uc},:) = Psth(Celltypei{uc},:) ;
    Psth_bytypes_toPsthi(Celltypei{uc}) = Psth_bytypei{uc} ;
    
    c = c+length(Celltypei{uc}) ;
end

%Psth = Psth_bytypes ; % OVERWRITE!!!!!

% % autocorr psth 
% for c = 1:NumCells ; % for each cell
%     PsthAc(c,:) = xcov(Psth(c,:)) ;
% end

% make psth and movClip 'Blocks

% define pnts you want to use for each block

% % automatically make some block pnts - SKIP THIS IF YOU JUST WANT TO MANUALLY DO IT
% % number of bins per psth and mov block 
% BlockBinTime = 1200 ; % (sec) aprox time span you want the blocks to be
% BlocksToMake = [2] ; % psth and movClip blocks you want to create and save, all blocks if empty
% BlockBinNumber = floor(BlockBinTime/diff(PsthTime(1:2))) ; % number of bins
% 
% if isempty(BlocksToMake) ; % if you want to make all possible blocks
%     BlocksToMake = [1:ceil(size(mov,3)/BlockBinNumber)] ;
% end
% 
% for block=1:length(BlocksToMake) ; % for each block you want to make
%     endPnt = min(length(Psth),BlocksToMake(block)*BlockBinNumber) ; % end point so last block can be smaller
%     BlockPnts{block} = [(BlocksToMake(block)-1)*BlockBinNumber+1:endPnt] ;
% end
% 
% % define blocks by psth-mov time pnts
% BlockPnts{1} = [1:75000] ;
% BlockPnts{2} = [76000:106000] ;
% BlockPnts{3} = [106001:106010] ;
% 
% % psth block
% for block=1:length(BlockPnts) ; % for each block 
%     PsthBlock{block} = Psth(:,BlockPnts{block}) ;
% end
% 
% % movie clip blocks (set of frames (ClipFrames in length) preceding mov time bin) in blocks    
% ClipFrames = 15 ; % number of movie frames 
% for block=1:length(BlockPnts) ; % for each block  
%     movClips{block} = zeros([size(mov,1),size(mov,2),length(BlockPnts{block}),ClipFrames]) ; % prep clip block
%     
%     for mc = 1:length(BlockPnts{block}) ; % for each frame
%         sf = max(BlockPnts{block}(mc)-ClipFrames+1,1) ; % mov start frame
%         sf2 = min(BlockPnts{block}(mc),ClipFrames) ; % movClips start frame
%         movClips{block}(:,:,mc,end-sf2+1:end) = mov(:,:,sf:BlockPnts{block}(mc)) ;
%     end
% end
% 
% % rearrange movClips so that its size(ClipFrame,Y,X,Frames) instead of (frames, X,Y, ClipFrame)
% for block=1:length(BlockPnts) ; % for each block  
%     movClips{block} = permute(movClips{block},[4,2,1,3]) ;
% end
% 
% % get rid of arrays for python
% for block=1:length(BlockPnts) ; % for each block
%     eval(['PsthBlock_',num2str(block),'= PsthBlock{block} ;']) ;
%     eval(['movClips_',num2str(block),'= int8(movClips{block}) ;']) ;
% end

% save mov and psth
save('/Users/jcafaro/Desktop/BWDb1.mat','mov','Psth','-v7.3') 

% % make STA (for sanity check)
% TestCelli = 12 ; % test cell
% maxLag = 20 ;
% Sta = nans(size(mov,1),size(mov,2),maxLag) ;
% for lag = 1:maxLag ; % for each time lag
%     for X = 1:size(mov,1) ; % proceed pix by pixel along X
%         for Y = 1:size(mov,2) ; % proceed pix by pixel along Y
%             Sta(X,Y,lag) = [Psth(TestCelli,lag:end),zeros(1,lag-1)] * squeeze(mov(X,Y,:)); % weighted sum in time 
%         end
%     end
% end

% % make LN prediction and assess corr coef on help out data
% blank = zeros(1,size(dataRun.stas.stas{1},4)) ;
% gen = zeros(NumCells,size(mov,3)) ; % prep
% for c = 1:NumCells ; % for each cell
%     disp(c)
%     for X = 1:size(mov,1) ;
%         for Y = 1:size(mov,2) ;
%             filt = [blank,fliplr(squeeze(dataRun.stas.stas{c}(X,Y,1,:))')] ; % filter
%             gen(c,:) = gen(c,:) + conv(squeeze(mov(X,Y,:)),filt,'same')' ;
%         end
%     end
%     
%     [Nlx(c,:),Nly(c,:)] = NonLinFilterFinder(gen(c,:),smooth(Psth(c,:),3)',range(gen(c,:))/200) ; % binned NL interpolated across missing data
%     LnEst(c,:) =  NonLinFilter(gen(c,:),Nlx(c,:),Nly(c,:)) ; % implement NL
%     
%     temp = corrcoef(smooth(Psth(c,:),3),LnEst(c,:)) ;
%     pearson(c) = temp(2) ;
% end
 
% LN prediction using significant stixels in STA for noise reduction
for c = 1:NumCells ; % for each cell
    dataRun = get_snls(dataRun, dataRun.cell_ids(c)) ;
    
    if sum(dataRun.stas.snls{c}.gen_signal).0 ;
        [Nlx2(c,:),Nly2(c,:)] = NonLinFilterFinder(dataRun.stas.snls{c}.gen_signal',full(dataRun.stas.snls{c}.spikes)',range(dataRun.stas.snls{c}.gen_signal)/200) ; % binned NL interpolated across missing data
        LnEst2(c,:) =  NonLinFilter(dataRun.stas.snls{c}.gen_signal,Nlx2(c,:),Nly2(c,:)) ; % implement NL

        temp = corrcoef(smooth(full(dataRun.stas.snls{c}.spikes),3),LnEst2(c,:)) ;
        pearsonLN(c) = temp(2) ;
    end
end

%% Repeating BW data
dataRunBwRep = load_data(DataBlock(DB).BwRepPath{Params.BwRepPathNum}) ;
dataRunBwRep = load_neurons(dataRunBwRep) ;

% get triggers
FramesPerTrig = 100 ; % frames per triger
FramesPerSec = 60 ; % frames per second
NumTrigsPerRep = round(Params.RepDuration*FramesPerSec/FramesPerTrig) ; % number of triggers per rep
RepTrigs = dataRunBwRep.triggers(1:NumTrigsPerRep:end) ;

% get movie
dataRunBwRep = load_java_movie(dataRunBwRep,DataBlock(DB).BwRepMoviePath{Params.BwRepPathNum}) ;
movTempBwRep = get_movie(DataBlock(DB).BwRepMoviePath{Params.BwRepPathNum}, dataRunBwRep.triggers, size(dataRunBwRep.stimulus.java_movie)) ;
movBwRep(:,:,:) = movTempBwRep(:,:,1,1:end-100) ; % get rid of rgb for achromatic stim and clip set frames after last trigger 

% make time bins 
TimeStep = 1/FramesPerSec ; % (sec)
trigdiff = diff(RepTrigs) ;
PsthTime = [0:TimeStep:diff(RepTrigs(1:2))] ;
  
% make psth
NumCells = length(dataRun.spikes) ; % for each cell in the master

for c = 1:NumCells ; % for each cell
    PsthBwRep{c} = nan(length(RepTrigs),length(PsthTime)-1) ;
    if sum(dataRunBwRep.cell_ids==dataRun.cell_ids(c))>0 ;
        celli = get_cell_indices(dataRunBwRep,dataRun.cell_ids(c)) ;
        for t=1:length(RepTrigs) ; % for each trigger
            TempBins = PsthTime+RepTrigs(t) ;
            PsthBwRep{c}(t,:) = histcounts(dataRunBwRep.spikes{celli},TempBins)./TimeStep ; % spikes per sec
        end
        PsthBwRep_mean(c,:) = mean(PsthBwRep{c}) ; % average across trials
    end
end

% % pearson coef (fidelidy) 
% for c = 1:NumCells ; % for each cell
%     for t=1:length(RepTrigs) ; % for each trigger
%         temp = corrcoef(PsthBwRep_mean(c,:),PsthBwRep{c}(t,:)) ;
%         PsthBwRep_r(c,t) = temp(2) ;
%     end
%     PsthBwRep_r_mean(c) = mean(PsthBwRep_r(c,:)) ; % avearge across trials
% end

% pearson coef (fidelidy) - estimate with averages 
NumResamples = 20 ;
for c = 1:NumCells ; % for each cell
    for t=1:NumResamples ; % for each trigger
        i1 = randperm(length(RepTrigs),length(RepTrigs)/2) ; %grab half the trials
        Psth_av1 = mean(PsthBwRep{c}(i1,:)) ;
        Psth_av2 = mean(PsthBwRep{c}(~ismember([1:length(RepTrigs)],i1),:)) ;
        temp = corrcoef(Psth_av1,Psth_av2) ;
        PsthBwRep_r(c,t) = temp(2) ;
    end
    PsthBwRep_r_mean(c) = mean(PsthBwRep_r(c,:)) ; % avearge across trials
end

% test LN model on BW repeat
dataRunBwRep.stas = dataRun.stas ;
dataRunBwRep.stas = rmfield(dataRunBwRep.stas,'snls') ;
for c = 1:NumCells ; % for each cell
    if sum(dataRunBwRep.cell_ids==dataRun.cell_ids(c))>0 ;
        dataRunBwRep = get_snls(dataRunBwRep, dataRun.cell_ids(c)) ;
        ci = get_cell_indices(dataRunBwRep, dataRun.cell_ids(c)) ;
        
        BwRepLnEst(c,:) =  NonLinFilter(dataRunBwRep.stas.snls{c}.gen_signal,Nlx2(c,:),Nly2(c,:)) ; % implement NL

        temp = corrcoef(dataRunBwRep.stas.snls{c} ,BwRepLnEst(c,:)) ;
        pearson2(c) = temp(2) ;
    end
end


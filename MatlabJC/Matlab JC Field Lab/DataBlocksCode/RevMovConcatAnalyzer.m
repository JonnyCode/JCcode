function ForIgor = RevMovConcatAnalyzer(DataBlock, DB, Params)

% analayze concatinated sets of response to movie, followed by movie played
% backwards

% JC 2018-02-07

% parameters
sampleRate = 20000 ;
PsthBinStepTime = .001 ; % (sec)
PsthBinSmoothTime = .020 ; % (sec)

NumMovs = length(DataBlock(DB).ConcatRevMoviePath) ;
NumReps = DataBlock(DB).ConcatRevMovieNumReps ;

PsthBinSmoothPnts = PsthBinSmoothTime*sampleRate ;

for m = 1:NumMovs ; % for each movie
    
    % load data
    dataRun = load_data(DataBlock(DB).ConcatRevMoviePath{m}) ;
    dataRun = load_neurons(dataRun) ;
    
    TrigsPerMovie = length(dataRun.triggers)/(2*NumReps(m)) ; % number of triggers per movie
    TrigsFwd = [1:TrigsPerMovie:length(dataRun.triggers)/2] ; % triggers of movie fwd
    TrigsRvs = [((length(dataRun.triggers)/2)+1):TrigsPerMovie:length(dataRun.triggers)] ; % triggers of movie in rvs
    
    PsthTime(m,:) = [0:PsthBinStepTime:diff(TrigsFwd(1:2))] ;
    
    NumCells = length(dataRun.spikes) ;
    for c = 1:NumCells ; % for each cell
        for t = 1:NumReps(m) ; % for each Rep
            spks = dataRun.spikes{c}-TrigsFwd(t) ;
            spkTrainFwd(t,:) = histcounts(spks,PsthTime(m,:)) ;
            
            spks = dataRun.spikes{c}-TrigsRvs(t) ;
            spksRvs = -(PsthTime(m,end) - spks) ;
            spkTrainRvs(t,:) = histcounts(spksRvs,PsthTime(m,:)) ;
        end
        PsthFwd(c,:) = mean(spkTrainFwd,1) ;
        PsthRvs(c,:) = mean(spkTrainRvs,1) ;
        
        PsthFwd_smooth(c,:)= smooth(PsthFwd(c,:),PsthBinSmoothPnts) ;
        PsthRvs_smooth(c,:) = smooth(PsthRvs(c,:),PsthBinSmoothPnts) ;
        
        r(c) = corr(PsthFwd_smooth(c,200:end-200)',PsthRvs_smooth(c,200:end-200)') ;
    end
end

    
    


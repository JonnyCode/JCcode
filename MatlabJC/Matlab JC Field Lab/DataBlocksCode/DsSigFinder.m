function ForIgor = DsSigFinder(DataBlock,DB,Params)

% JC 4/30/2018
% this function will find cells with significant response distributions to
% drifiting grating

% NEVER VALIDATED!!!!! - seemed to select too many DS cells


TrialTrigInterval = 10 ; % sec

% find ds cells with drifting gratings data

% load data 

if isfield(Params,'TimeBounds') ; % if you are selecting only part of dataRun
    dataRun = load_data(DataBlock(DB).DsConcatPath) ;
    dataRun = load_neurons(dataRun) ;
    dataRun.triggers  = dataRun.triggers(dataRun.triggers>=Params.TimeBounds(1) & dataRun.triggers<=Params.TimeBounds(2)) ; 
elseif isfield(Params,'DsPathNum') ;
    dataRun = load_data(DataBlock(DB).DsPath{Params.DsPathNum}) ; 
    dataRun = load_neurons(dataRun) ;
else 
    dataRun = load_data(DataBlock(DB).DsPath) ;
    dataRun = load_neurons(dataRun) ;
end

dataRun = load_ei(dataRun, 'all') ;

% stim path
if isfield(Params,'DsPathNum') ;
    slashi = strfind(DataBlock(DB).DsPath{Params.DsPathNum},'/') ; % find the /
    dataRun.names.stimulus_path = [DataBlock(DB).DsPath{Params.DsPathNum}(1:slashi(end-1)),'stimuli/s',DataBlock(DB).DsPath{Params.DsPathNum}(end-1:end),'.txt'] ;
else
    slashi = strfind(DataBlock(DB).DsPath,'/') ; % find the /
    dataRun.names.stimulus_path = [DataBlock(DB).DsPath(1:slashi(end-1)),'stimuli/s',DataBlock(DB).DsPath(end-1:end),'.txt'] ;  
end
    
dataRun = load_stim(dataRun,'user_defined_trigger_interval', TrialTrigInterval) ;

%parameters
NumCells = length(dataRun.spikes) ;
NumDirDg = length(dataRun.stimulus.params.DIRECTION) ; % number of directions
NumTempPeriods = length(dataRun.stimulus.params.TEMPORAL_PERIOD) ; % number of temporal periods

% stim trials
for t=1:length(dataRun.stimulus.trials) ; % for each DG trial
    DgParams(t,1) = dataRun.stimulus.trials(t).DIRECTION ; % direction
    DgParams(t,2) = dataRun.stimulus.trials(t).TEMPORAL_PERIOD; % temporal period
end 

% Dg spike count
for cells=1:NumCells ; % for each cell
    for st = 1:NumDirDg ; % for each direction 
        for tmp = 1:NumTempPeriods ; % for each temp period
            prm = [dataRun.stimulus.params.DIRECTION(st),dataRun.stimulus.params.TEMPORAL_PERIOD(tmp)] ; % params set
            ti = find((DgParams(:,1)==prm(1)).*(DgParams(:,2)==prm(2))==1) ; % index of triggers for that stim
            for t=1:length(ti) ; % for each trial
                spk = dataRun.spikes{cells}-dataRun.stimulus.triggers(ti(t)) ;
                spk = spk(spk>=0 & spk<=TrialTrigInterval) ;

                spkCount(cells,st,tmp,t) = length(spk) ;
            end
            
        end
        spkCount_mean(cells,st) = mean(spkCount(cells,st,:)) ; % average across trials and tmp
    end
end

% grating vector sum X,Y
for cells=1:NumCells ; % for each cell

    for st = 1:NumDirDg ; % for each stimulus 
        VectTempX(st) = cosd(dataRun.stimulus.params.DIRECTION(st))*spkCount_mean(cells,st) ; 
        VectTempY(st) = sind(dataRun.stimulus.params.DIRECTION(st))*spkCount_mean(cells,st) ; 
    end
    VectSumX(cells) = sum(VectTempX) ;
    VectSumY(cells) = sum(VectTempY) ;
end 

% projection onto each trial
for cells=1:NumCells ; % for each cell
    TempDot = [] ; % empty vect
    for tmp = 1:NumTempPeriods ; % for each temp period
        for t=1:length(ti) ; % for each trial
            VectTempX = sum(cosd(dataRun.stimulus.params.DIRECTION).*spkCount(cells,:,tmp,t)) ; 
            VectTempY = sum(sind(dataRun.stimulus.params.DIRECTION).*spkCount(cells,:,tmp,t)) ;
            TempDot = [TempDot,[VectTempX,VectTempY]*[VectSumX(cells);VectSumY(cells)]] ;
        end
    end
    [Sig(cells),p(cells)] = ttest(TempDot) ; % not zero
end

    
            
           
        






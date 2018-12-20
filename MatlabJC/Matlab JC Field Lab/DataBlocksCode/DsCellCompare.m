function ForIgor = KoDsConcatAnalysis(DataBlock, DB, Params)

% this function will compare DS cell DG responses from a concatinated datablock
% JC 2017-12-04

DsPathNums = Params.DsPathNums ; % path numbers to compare
ConcatPathNum = 1 ;

% parameters
Color_list = {'k','r','b','g','c','y','k','r','b','g','c','y','k','r','b','g','c','y','k','r','b','g','c','y'} ; % order of colors for each 

saveDsIdsPath = ['/Users/jcafaro/Documents/AnalysisFigures/KoCompare/DsIdsDb',num2str(DB),'ConcatPathNum',num2str(ConcatPathNum)] ;
 % stim path
TrialTrigInterval = 10 ;
slashi = strfind(DataBlock(DB).DsPath{DsPathNum},'/') ; % find the /
dataRunDg.triggers = dataRun.triggers(dataRun.triggers<=DgStimTime) ;
dataRunDg.names.stimulus_path = [DataBlock(DB).DsPath{Params.DsPathNum}(1:slashi(end-1)),'stimuli/s',DataBlock(DB).DsPath{DsPathNum}(end-1:end),'.txt'] ;
dataRunDg = load_stim(dataRunDg,'user_defined_trigger_interval', TrialTrigInterval) ;

% psth
%psthTimeDg = [0:StimFrameRate:TrialTrigInterval] ; 
psthTimeDg = [0:psthBinTime:TrialTrigInterval] ;
for cells=1:length(dataRun.spikes) ; % for each cell
    for st = 1:length(dataRunDg.stimulus.combinations) ; % for each stimulus 
        ti = find(dataRunDg.stimulus.trial_list==st) ; % index of triggers for that stim
        for t=1:length(ti) ; % for each trial
            spk = dataRun.spikes{cells}-dataRunDg.stimulus.triggers(ti(t)) ;
            spk = spk(spk>=0 & spk<=TrialTrigInterval) ;
            %psthDg{cells}{st}(t,:) = histc(spk,psthTimeDg) ;
            
            spikeTrain = zeros(1,length(psthTimeDg)) ; % empty spike train
            spkPnts = 1+floor(spk/psthBinTime) ; % spike points 
            spikeTrain(spkPnts) = 1 ; % spike train
            psthDg{cells}{st}(t,:) = smooth(spikeTrain,psthSmoothPnts) ; % psth
        end
        psthDg_mean{cells}(st,:)=mean(psthDg{cells}{st},1) ;
        psthDg_var{cells}(st,:)=var(psthDg{cells}{st},[],1) ;
    end
end

% grating vector sum
for cells=1:length(dataRun.spikes) ; % for each cell
    VectTemp = nans(length(dataRunDg.stimulus.combinations),2) ; % prep mat
    for st = 1:length(dataRunDg.stimulus.combinations) ; % for each stimulus 
        VectTemp(st,:) = [dataRunDg.stimulus.combinations(st).DIRECTION,...
            sum(psthDg_mean{cells}(st,ceil(StimOnOffset/diff(psthTimeDg(1:2))):ceil(DgOffTime/diff(psthTimeDg(1:2)))))] ; % prep vector (dir, spike rate)
    end
    VectSumDg(cells,:) = PolarVectorAddition(VectTemp) ;
end      

% grating tuning curves divided across speeds
for st = 1:length(dataRunDg.stimulus.combinations) ; % for each stimulus 
    di = find(dataRunDg.stimulus.params.DIRECTION == dataRunDg.stimulus.combinations(st).DIRECTION) ;
    ti = find(dataRunDg.stimulus.params.TEMPORAL_PERIOD == dataRunDg.stimulus.combinations(st).TEMPORAL_PERIOD) ;
    for cells=1:length(dataRun.spikes) ; % for each cell
        DgTuningCurve{ti}(cells,di) = sum(psthDg_mean{cells}(st,ceil(StimOnOffset/diff(psthTimeDg(1:2))):ceil(DgOffTime/diff(psthTimeDg(1:2))))) ; %  ; % NOT FINISHED
    end
end
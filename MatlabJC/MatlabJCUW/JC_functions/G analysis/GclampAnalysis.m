% this script will analyze current clamp data from cells presented 
% repeating blocks of mulitple random light noise stimuli, g clamp, and posiible
% recorded cell attached

% JC 6/23/08


%% input ca, exc, and inh epochs and pre and stim points 

CellName_str = '070908Bc1' ;

Spike_epochs_str = {'[25:44]'} ;

Iclamp_epochs_str = {'[88:102,229:238]'} ;

Gclamp_epochs_str = {'[159:168,194:203]'} ;

PrePnts = 10000 ;   % number of points before stimuli starts
StmPnts = 60000 ;   % number of points during stimuli
STA_pnts = 3000 ;   % number of points used for prespike wave forms


%% get light stimuli and data for all epochs
[fp, error] = ITCInitializeAnalysis(500000, ['~/Data/Primate/',CellName_str]);

Spike_epochs = str2num(Spike_epochs_str{1}) ;
for a = 1:length(Spike_epochs) ; % for each spike epoch
    [Spike_stm(a,:), error] = ITCReadEpochStm(Spike_epochs(a), 0, fp);  % get the light stimulus
    [Spike_data(a,:), error] = ITCReadEpoch(Spike_epochs(a), 0, fp);    % get cell attached data
end

Iclamp_epochs = str2num(Iclamp_epochs_str{1}) ;
for a = 1:length(Iclamp_epochs) ; % for each current clamp epoch
    [Iclamp_stm(a,:), error] = ITCReadEpochStm(Iclamp_epochs(a), 0, fp);      % get the light stimulus
    [Iclamp_data(a,:), error] = ITCReadEpoch(Iclamp_epochs(a), 0, fp);        % get current clamp data
    [SamplingInterval_IC(a), error] = ITCGetSamplingInterval(Iclamp_epochs(a), fp);        % get sample interval    
    SamplingInterval_IC(a) = SamplingInterval_IC(a) * 1e-6;                                   % change us to seconds
end

Gclamp_epochs = str2num(Gclamp_epochs_str{1}) ;
for a = 1:length(Gclamp_epochs) ; % for each conductance clamp epoch
    [Gclamp_exc(a,:), Gclamp_inh(a,:), error] = ITCReadEpochStmGClamp(Gclamp_epochs(a), 0, fp);    % get the conductances used
    [Gclamp_data(a,:), error] = ITCReadEpoch(Gclamp_epochs(a), 0, fp);                  % get gclamp data
    [SamplingInterval_GC(a), error] = ITCGetSamplingInterval(Gclamp_epochs(a), fp);        % get sample interval    
    SamplingInterval_GC(a) = SamplingInterval_GC(a) * 1e-6;                                   % change us to seconds
end

% check that stimuli for ca and Iclamp are all in the same order
if ~isequal(Spike_stm,Iclamp_stm) ; % if all the stim used are in the correct order 
    print('data out of order?')
end

% subtract offsets from the cell attached data only 
Spike_data = Spike_data - repmat(mean(Spike_data(:,1:PrePnts),2),1,size(Spike_data,2)) ; % subtract mean of prepoints

% detect spikes in cell attached data, current clamp data, and coductance
% clamp data

SpikePnts_CA = SpikeDetection(Spike_data,10,10000) ; % data,threshold,samplerate
SpikePnts_IC = SpikeDetection_WC(Iclamp_data,-20,1/SamplingInterval_IC(1)) ; % data, threshold, samplerate
SpikePnts_GC = SpikeDetection_WC(Gclamp_data,-20,1/SamplingInterval_GC(1)) ; % data, theshold, samplerate


%% get pre iclamp spike light and voltage

for a = 1:length(SpikePnts_IC) ;                                                                    % for each spike epoch ...
    STA_spikesIC{a} = find(SpikePnts_IC{a}>PrePnts+STA_pnts & SpikePnts_IC{a}<PrePnts+StmPnts) ;        % indicies of SpikePnts vector during the stimuli and far enough out to get the prespike wave form 
    PreSpikeIC_stm{a} = nans(length(STA_spikesIC{a}),STA_pnts) ;                                        % prep the matrix for speed
    PreSpikeIC_volt{a} = nans(length(STA_spikesIC{a}),STA_pnts) ;
    
    for b = 1:length(STA_spikesIC{a}) ;                                                       % for each spike that can be used to note a prespike waveform...                                      
            PreSpikeIC_stm{a}(b,:) = Iclamp_stm(a,SpikePnts_IC{a}(STA_spikesIC{a}(b))-STA_pnts+1:SpikePnts_IC{a}(STA_spikesIC{a}(b))) ;
            PreSpikeIC_volt{a}(b,:) = Iclamp_data(a,SpikePnts_IC{a}(STA_spikesIC{a}(b))-STA_pnts+1:SpikePnts_IC{a}(STA_spikesIC{a}(b))) ;                        
    end % spike loop
    
    PreSpikeIC_stm_mean(a,:) = mean(PreSpikeIC_stm{a}) ; % mean prespike waveform for each epoch
    PreSpikeIC_volt_mean(a,:) = mean(PreSpikeIC_volt{a}) ;
    
    if a == 1 ; % if this is the first epoch
        PreSpikeIC_stm_all = PreSpikeIC_stm{a} ;
        PreSpikeIC_volt_all = PreSpikeIC_volt{a} ;
    else
    PreSpikeIC_stm_all = [PreSpikeIC_stm_all ; PreSpikeIC_stm{a}] ; % concatinate prespike waveforms from all epochs
    PreSpikeIC_volt_all = [PreSpikeIC_volt_all ; PreSpikeIC_volt{a}] ;
    end
end % end epoch loop

% get average pre spike waveforms over all epochs   
STAIC_stm = mean(PreSpikeIC_stm_all) ;
STAIC_volt = mean(PreSpikeIC_volt_all) ;

%% get pre Gclamp spike light(from iclamp data) and voltage

for a = 1:length(SpikePnts_GC) ;                                                                    % for each spike epoch ...
    SpikePnts_GCrs{a} = round(SpikePnts_GC{a} * (SamplingInterval_GC(a)/SamplingInterval_IC(a))) ;  % convert gc pnts into ic points
    STA_spikesGCrs{a} = find(SpikePnts_GCrs{a}>PrePnts+STA_pnts & SpikePnts_GCrs{a}<PrePnts+StmPnts) ;        % indicies of SpikePnts vector during the stimuli and far enough out to get the prespike wave form 
    PreSpikeGCrs_stm{a} = nans(length(STA_spikesGCrs{a}),STA_pnts) ;                                        % prep the matrix for speed
    
    for b = 1:length(STA_spikesGCrs{a}) ;                                                       % for each spike that can be used to note a prespike waveform...                                      
            PreSpikeGCrs_stm{a}(b,:) = Iclamp_stm(a,SpikePnts_GCrs{a}(STA_spikesGCrs{a}(b))-STA_pnts+1:SpikePnts_GCrs{a}(STA_spikesGCrs{a}(b))) ;
                              
    end % spike loop
    
    PreSpikeGCrs_stm_mean(a,:) = mean(PreSpikeGCrs_stm{a}) ; % mean prespike waveform for each epoch
    
    if a == 1 ; % if this is the first epoch
        PreSpikeGCrs_stm_all = PreSpikeGCrs_stm{a} ;
    else
    PreSpikeGCrs_stm_all = [PreSpikeGCrs_stm_all ; PreSpikeGCrs_stm{a}] ; % concatinate prespike waveforms from all epochs
    end
end % end epoch loop

% get average pre spike waveforms over all epochs   
STAGCrs_stm = mean(PreSpikeGCrs_stm_all) ;

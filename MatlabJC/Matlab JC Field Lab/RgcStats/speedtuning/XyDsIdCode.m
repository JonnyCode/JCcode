
DB = 21 ; % 19,20,21
[DataBlock,Params] = DataBlocks_NaturalStim ;

DsPathNum = 2 ;
TrialTrigInterval = 10 ;


% load data
dataRun = load_data(DataBlock(DB).DsPath{DsPathNum}) ;
dataRun = load_neurons(dataRun) ;
dataRun.piece.array_id = 1551 ; % CAUTION

if isfield(Params,'DsPathNum') ;
    slashi = strfind(DataBlock(DB).DsPath{DsPathNum},'/') ; % find the /
    dataRun.names.stimulus_path = [DataBlock(DB).DsPath{DsPathNum}(1:slashi(end-1)),'stimuli/s',DataBlock(DB).DsPath{DsPathNum}(end-1:end),'.txt'] ;
else
    slashi = strfind(DataBlock(DB).DsPath,'/') ; % find the /
    dataRun.names.stimulus_path = [DataBlock(DB).DsPath(1:slashi(end-1)),'stimuli/s',DataBlock(DB).DsPath(end-1:end),'.txt'] ;  
end
    
dataRun = load_stim(dataRun,'user_defined_trigger_interval', TrialTrigInterval) ;



%% separate ds vs non-ds cells

[NumSpikesCell, StimComb] = get_spikescellstim(dataRun,dataRun.cell_ids,0);

ds_struct = dscellanalysisWithFigs(NumSpikesCell, StimComb); % plot vector sums for each speed

params_idx = [2 3]; % which speeds to use for classification

[ds_id, nonds_id] = classify_ds(dataRun, ds_struct, params_idx);



%% separate ON from ON-OFF by speed

[id_sub, idx_sub] = classify_onoff_from_speed(datarun, ds_id, pc1, pc2, manual)
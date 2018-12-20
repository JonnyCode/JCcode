function ForIgor = DuplicateSearch(DataBlock, DB, Params) 

% look for duplicates

% JC 8/19/16

opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1,'load_sta',1) ;
dataRun = load_data(DataBlock(DB).BW,opt) ;

Duplicate_cell_ids = DuplicateFinder(dataRunMaster,dataRunSlave,varargin)
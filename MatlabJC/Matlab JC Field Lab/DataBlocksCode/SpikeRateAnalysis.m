function ForIgor = SpikeRateAnalysis(DataBlock, DB, Params)

% This function will examine the spontaneous spike rate and change during bwn 

% JC 6/19/15

% parameters
Color_list = {'c','k','r','g','y','b'} ; % order of colors for each 

% load data
dataRun = load_data(DataBlock(DB).SpontPath{1}) ;
dataRun = load_neurons(dataRun) ;
dataRun = load_params(dataRun) ;

% params
NumCells = length(dataRun.spikes) ;

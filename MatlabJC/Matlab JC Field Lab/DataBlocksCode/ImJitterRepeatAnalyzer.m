function ForIgor = ImJitterRepeatAnalyzer(DataBlock, DB, Params)

% analyze ImJitterRepeat data

%JC 2018-11-30

% data set select
RunAsScript = true ;
if RunAsScript ;
    DB = 21 ; % 19,20,21
    [DataBlock,Params] = DataBlocks_NaturalStim ;
end

DsPathNum = 1 ;

% flags and parameters 
Color_list = {'c','r','b','g','k','y','c','r','b','g','k','y','c','r','b','g','k','y',...
    'c','r','b','g','k','y','c','r','b','g','k','y','c','r','b','g','k','y'} ; % order of colors for each 
FramesPerTrig = 100 ; % number of frames between each trigger

% save and loads paths
saveDsIdsPath = ['/Users/jcafaro/Documents/AnalysisFigures/NatStimDs/ImJitterRepeatAnalyzer/DsIdsDb',num2str(DB),'DsPathNum',num2str(DsPathNum)] ;

% load data
dataRun = load_data(DataBlock(DB).ImJitterRepeats) ;
dataRun = load_neurons(dataRun) ;
dataRun.piece.array_id = 1551 ; % CAUTION - assumes dense array
dataRun = load_ei(dataRun, 'all') ;

% identify DS cells
Params.OutlierFlag = true ; 
Params.DsPathNum = DsPathNum ;
Params.DataBlockToMap = DataBlock(DB).ImJitterRepeats ;

if UseImportDsIdsPath ; % if you want the DS ids selected elsewhere
    load(ImportDsIdsPath) ;
else
    try load(saveDsIdsPath) ; % if they are ds ids already saved
    catch % if not find them
        ForIgor = DsCellFinder(DataBlock, DB, Params) ;
        save(saveDsIdsPath, 'ForIgor')
    end
end

cell_id = ForIgor.master_id{2} ; % id of repeated data

for DsType=1:length(cell_id) ; % cell indicies of DsTypes
    cell_i{DsType} = get_cell_indices(dataRun, cell_id{DsType}) ;
end

% load stimulus
slashi = strfind(DataBlock(DB).ImJitterRepeats,'/') ; % find the /
dashi = strfind(DataBlock(DB).ImJitterRepeats,'-') ; % find the -
StimPath = [DataBlock(DB).ImJitterRepeats(1:slashi(end-1)),'stimuli/s0',DataBlock(DB).ImJitterRepeats(end)] ;
load(StimPath) ;

% number of triggs for each repeat
NumTriggsPerRep = stimulus.num_frames/FramesPerTrig ;

Trigs = dataRun.triggers(1:NumTriggsPerRep:end) ; % triggers at begining of each repeat

%% figure
RunId = num2str(randi(100000)) ; % random number to save figs and mat files for this code run

for DsType=1:length(cell_i) ; % cell indicies of DsTypes
    for c=1:length(cell_i{DsType}) ;
        figure(1) % raster of ds cells
        clf
        spikes_by_trials = get_raster(dataRun.spikes{cell_i{DsType}(c)}, Trigs,...
            'axis_range',[50,70,0,stimulus.num_repeats+1],'dots_flag',false) ;
        pause
    end
end

saveas(gcf,[saveFigPath,'Raster',num2str(DB),'Fig_',RunId])
print(gcf, '-dpdf',[saveFigPath,'Raster',num2str(DB),'Fig_',RunId])
 

%%


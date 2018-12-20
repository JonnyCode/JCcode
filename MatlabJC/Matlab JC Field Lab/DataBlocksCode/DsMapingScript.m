% script to get ds cell ids and map onto different data sets for ImjitterV4 code

MapSetNum = 2 ;

% import datablocks
[DataBlock,Params] = DataBlocks_NaturalStim ;

% ds id mapped sets
MapSet(1).DB = 19 ;
MapSet(1).DsPathNum = 1 ;
MapSet(1).DataBlockToMap{1} = DataBlock(19).ImJitter{1} ; 
MapSet(1).DataBlockToMap{2} = DataBlock(19).ImJitter{2} ;

MapSet(2).DB = 20 ;
MapSet(2).DsPathNum = 1 ;
MapSet(2).DataBlockToMap{1} = DataBlock(20).ImJitter{1} ; 
MapSet(2).DataBlockToMap{2} = DataBlock(20).ImJitter{2} ;

MapSet(3).DB = 21 ;
MapSet(3).DsPathNum = 1 ;
MapSet(3).DataBlockToMap{1} = DataBlock(21).ImJitter{1} ; 
MapSet(3).DataBlockToMap{2} = DataBlock(21).ImJitter{2} ;

MapSet(4).DB = 21 ;
MapSet(4).DsPathNum = 1 ;
MapSet(4).DataBlockToMap{1} = DataBlock(21).ImJitter{1} ; 

% ds id path
saveDsIdsPath = ['/Users/jcafaro/Documents/AnalysisFigures/',...
    'NatStimDs/ImJitterAnalysisV4DsSelection/DsIdsDb',...
    num2str(MapSet(MapSetNum).DB),'DsPathNum',num2str(MapSetNum)] ;

% find ds cells
Params.OutlierFlag = true ; 
Params.DsPathNum = MapSet(MapSetNum).DsPathNum ;

try load(saveDsIdsPath) ; % if they are ds ids already saved
catch % if not find them
    ForIgor = DsCellFinder(DataBlock, MapSet(MapSetNum).DB, Params) ;
    save(saveDsIdsPath, 'ForIgor')
end

% mapping
ds_id = ForIgor.ds_id ; 
CommonMappedIds = ds_id ; % initialize

dataRun = load_data(DataBlock(MapSet(MapSetNum).DB).DsPath{MapSet(MapSetNum).DsPathNum}) ; 
dataRun = load_neurons(dataRun) ;
dataRun = load_ei(dataRun, 'all') ;

for db = 1:length(MapSet(MapSetNum).DataBlockToMap) ; % for each block to map

    dataRunMaster = load_data(MapSet(MapSetNum).DataBlockToMap{db}) ;
    dataRunMaster = load_neurons(dataRunMaster) ;
    dataRunMaster = load_ei(dataRunMaster, 'all') ;
    
    cell_ids{db} = dataRunMaster.cell_ids ;

    % map using electrical images
    cell_list_map = map_ei(dataRunMaster, dataRun) ; % ei map
    
    % change array to mat
    cell_list_map_mat{db} = nans(1,length(cell_list_map)) ; % prep mat
    for c=1:length(cell_list_map) ;
        if ~isempty(cell_list_map{c}) ;
            cell_list_map_mat{db}(c) = cell_list_map{c} ;
        end
    end
    
    % cells ids in slave for each UniqueCellType set in master dat
    for subi=1:2 ; % for each subset
        for drSelect=1:4 ; 
            [tempIds,mi] = intersect(cell_list_map_mat{db},ds_id{subi}{drSelect}) ; % ids of ds cells in dg set
            mapped_ds_ids{db}{subi}{drSelect} = tempIds ;
            MasterIds{db}{subi}{drSelect} = cell_ids{db}(mi) ; % ids of ds cells in master (jitter stim data)
            CommonMappedIds{subi}{drSelect} = intersect(CommonMappedIds{subi}{drSelect},mapped_ds_ids{db}{subi}{drSelect}) ; % dg set ids
        end
    end  
end

% cell ids in each set that are mapped in all sets
for db = 1:length(MapSet(MapSetNum).DataBlockToMap) ; % for each block to map
    for subi=1:2 ; % for each subset
        for drSelect=1:4 ; 
            [tempIds,mi] = intersect(cell_list_map_mat{db},CommonMappedIds{subi}{drSelect}) ; % ids of ds cells in dg set
            CommonMasterIds{db}{subi}{drSelect} = cell_ids{db}(mi) ; % ids of ds cells in master
        end
    end 
end



ImportDsIdsPath = ['/Users/jcafaro/Documents/AnalysisFigures/',...
    'NatStimDs/ImJitterAnalysisV4DsSelection/ImportDsIdsDb',...
    'MapSet',num2str(ImPathNum)] ;

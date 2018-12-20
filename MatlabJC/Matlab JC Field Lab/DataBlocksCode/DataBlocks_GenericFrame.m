% generic structure of DataBlocks script

%% functions to perform
perform.AnalysisX = false ;
perform.AnalysisY = false ;
perform.AnalysisZ = true ;

%% perform functions on this data block set
DBset_id = 1 ;

% data block sets
DBset{1} = [1] ; % temp analysis
DBset{2} = [1:3] ; % Cells with...

%% Data Blocks defined
% Root dircetory for analysis files
%RootAnalysisDir = '/Volumes/Berlioz/Analysis/' ; % Berlioz
RootAnalysisDir = '/Volumes/lab/Experiments/Array/Analysis/' ; % Brahms server
RootMovieDir = '/Volumes/lab/acquisition/movie-xml/' ; % Brahms server xml movie root

% VIP-cre mouse - bw stim
% Greg thought this one worked
DataBlock(1).PreKo.DataPath = [RootAnalysisDir,'2013-12-22-2/data000-3200-5000s-gdf/data000-3200-5000s'] ; % control
DataBlock(1).PostKo.DataPath = [RootAnalysisDir,'2013-12-22-2/data004/data004'] ; % +psem

% same retina as db 1 analyzed diffrently 
DataBlock(2).PreKo.DataPath = [RootAnalysisDir,'2013-12-22-2/data000-3200-5000s-gdf/data000-3200-5000s'] ; % control+10
DataBlock(2).PostKo.DataPath = [RootAnalysisDir,'2013-12-22-2/data004-data000-3200-5000s-gdf-map/data004-data000-3200-5000s-gdf-map'] ; % +psem
DataBlock(2).Wash.DataPath = [RootAnalysisDir,'2013-12-22-2/data006/data006'] ; % wash
DataBlock(2).mapEiFlag = false ; % don't map by ei

% Params defaults
Params.empty = [] ;
ForIgor = [] ;

%% run functions

if perform.AnalysisX ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = AnalysisX(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.AnalysisY ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = AnalysisY(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

%% population analysis



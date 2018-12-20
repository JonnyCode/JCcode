% DataBlocks for classification experiments

%% functions to perform
perform.SquareAnalysis = false ;

perform.DuplicateSearch = false ;

perform.SquareAnalysisConcat = false ;

perform.SquareAnalysisConcatWithBwMapping = true ;

%% perform functions on this data block set
DBset_id = 1 ;

% data block sets
DBset{1} = [2] ; % temp analysis
DBset{2} = [1:3] ; % Cells with...

%% Data Blocks defined
% Root dircetory for analysis files
%RootAnalysisDir = '/Volumes/Berlioz/Analysis/' ; % Berlioz
RootAnalysisDir = '/Volumes/lab/Experiments/Array/Analysis/' ; % Brahms server
RootMovieDir = '/Volumes/lab/acquisition/movie-xml/' ; % Brahms server xml movie root

% Thy1-egfp/yfp mouse on sparse array
% using photons for the first time
DataBlock(1).BwPath = [RootAnalysisDir,'2016-08-04-0/data000/data000'] ; % control (50% contrast higher mean)
DataBlock(1).BwPath = [RootAnalysisDir,'2016-08-04-0/data004/data004'] ; % control (100% contrast lower mean) 
DataBlock(1).MovingSquare = [RootAnalysisDir,'2016-08-04-0/data001/data001'] ; % 30 width, 10 shift
%DataBlock(1).MovingSquare = [RootAnalysisDir,'2016-08-04-0/data005/data005'] ; % 100 width, 25 shift
DataBlock(1).FfPulse = [RootAnalysisDir,'2016-08-04-0/data002/data002'] ; %
DataBlock(1).DsPath = [RootAnalysisDir,'2016-08-04-0/data003/data003'] ; % (2 TP, 8 directions)

% c57 mouse on dense array
% using photons (BW has oscillations)
DataBlock(2).BwPath{1} = [RootAnalysisDir,'2016-08-15-0/data000/data000'] ; % control 
DataBlock(2).MovingSquare = [RootAnalysisDir,'2016-08-15-0/data001/data001'] ; % 30 width, 10 shift
DataBlock(2).MovingSquare = [RootAnalysisDir,'2016-08-15-0/data002/data002'] ; % 100 width, 25 shift
DataBlock(2).FfPulse = [RootAnalysisDir,'2016-08-15-0/data003/data003'] ; %
DataBlock(2).DsPath = [RootAnalysisDir,'2016-08-15-0/data004/data004'] ; % (2 TP, 8 directions)
DataBlock(2).DataConcat = [RootAnalysisDir,'2016-08-15-0/data001-2-3-4/data001-2-3-4'] ; % concatinated data

%% Params defaults
Params.empty = [] ;
ForIgor = [] ;

%% run functions

if perform.SquareAnalysis ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = SquareAnalysis(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.SquareAnalysisConcat ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = SquareAnalysisConcat(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end


if perform.SquareAnalysisConcatWithBwMapping ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = SquareAnalysisConcatWithBwMapping(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end


%% population analysis




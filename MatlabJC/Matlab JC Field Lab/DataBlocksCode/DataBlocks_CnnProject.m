% DataBlocks for CNN project

%% functions to perform
perform.CnnPythonPrep = true ; 

%% perform functions on this data block set
DBset_id = 1 ;

% data block sets
DBset{1} = [1] ; % 

%% Root dircetory for analysis files
%RootAnalysisDir = '/Volumes/Berlioz/Analysis/' ; % Berlioz
RootAnalysisDir = '/Volumes/lab/Experiments/Array/Analysis/' ; % Brahms server
RootMovieDir = '/Volumes/lab/acquisition/movie-xml/' ; % Brahms server xml movie root
RootImageDir = '/Volumes/lab/Experiments/Array/Images/' ; 
RootSavePath = '/Users/jcafaro/Documents/AnalysisFigures/' ;
RootNatStimMoviesDir = '/Volumes/lab/Documents/Movies/' ; % natural stim moves

%% data blocks

% rat (NDF 0 before data005)
DataBlock(1).BwPath{1} = [RootAnalysisDir,'2017-01-16-0/data006_KR/data006_KR'] ; % 15-1
DataBlock(1).BwMoviePath{1} = [RootMovieDir,'BW-15-1-0.48-11111-53x40-60.35.xml'] ; % BW movie path
DataBlock(1).BwRepPath{1} = [RootAnalysisDir,'2017-01-16-0/data005-map_KR/data005-map_KR'] ; % 15-1 (NDF 0)
DataBlock(1).BwRepMoviePath{1} = [RootMovieDir,'BW-15-1-0.48-11111-53x40-60.35.xml'] ; % BW rep movie path
DataBlock(1).DsPath{1} = [RootAnalysisDir,'2017-01-16-0/data007-map/data007-map'] ; % Drifting grating
DataBlock(1).MovieRepPath{1} = [RootAnalysisDir,'2017-01-16-0/data009-map_KR/data009-map_KR'] ; % cat movie
DataBlock(1).NatMoviePath{1} = [RootNatStimMoviesDir,'CatCam/cat_mean117_sd62_0to255.mat'] ; % cat movie path
DataBlock(1).MovieStixWidth = 2 ; % movie stixel width
DataBlock(1).MovieFrameInterval = 1 ;% movie frame interval
DataBlock(1).MovieRepFrameNum = 300 ; % movie frame number

% Params
Params.BwPathNum = 1 ;
Params.BwRepPathNum = 1 ;
Params.displayFrameRate = 60.35 ;
Params.movieRefresh = 1 ;
Params.Movie_lag_time = 0 ;

% Params defaults
Params.empty = [] ;
ForIgor = [] ;

% run functions

if perform.CnnPythonPrep ;
    for a=1:length(DBset(DBset_id)) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = CnnPythonPrep(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

%% export to hdf5 file

% check that field names are not too long for Igor (must be less than 32 characters)
ForIgorFields = fieldnames(ForIgor) ;
for a=1:length(ForIgorFields) ; % for every field
    fieldLength(a) = length(ForIgorFields{a}) ;
end

% repository status
temp = getGitInfo ; 
RepVer = temp.hash ; % repository version 
[tempS,tempR] = system('git status') ;
if length(tempR)==62 ;
    ForIgor.RepStat = ['GitHub up to date ',RepVer] ; % if no unsynced files
else
    ForIgor.RepStat = ['GitHub not up to date ',RepVer] ;
end

if max(fieldLength)<32 ;
%     cd Z:/cafaro' Documents'/Analysis/ForIgorHdfs/   
%     exportStructToHDF5(ForIgor,'CNGfig.h5','/')
else disp('name too long for igor')
end


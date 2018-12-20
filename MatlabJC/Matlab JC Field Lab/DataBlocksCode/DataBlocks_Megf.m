% script to analyze CNG-/- +/- rescue, wt

% 3/17/2015 JC

% functions to perform
perform.MbAnalyzer = false ; 

perform.MBAnalyzerExtra = true ;

perform.DsCellBwAnalyzer = false ; 

% perform functions on this data block set
DBset_id = 1 ;

% data block sets
DBset{1} = [1:2] ; % 
DBset{2} = [1,2] ; % Het KO Megf10
DBset{3} = [3,4] ; % Hom KO Megf10
DBset{4} = [5] ; % wt mouse

% Root dircetory for analysis files
%RootAnalysisDir = '/Volumes/Berlioz/Analysis/' ; % Berlioz
RootAnalysisDir = '/Volumes/lab/Experiments/Array/Analysis/' ; % Brahms server
RootMovieDir = '/Volumes/lab/acquisition/movie-xml/' ; % Brahms server xml movie root
RootAnalysisSaveDir = '/Volumes/lab/Experiments/Array/Shared/JC/AnalysisMatFiles/' ; % local directory

% Megf10 homozygous KO
% includes image of stb
DataBlock(1).MbPath{1} = [RootAnalysisDir,'2015-09-01-0_RS/data001/data001'] ; % high contrast increments
DataBlock(1).MbPath{2} = [RootAnalysisDir,'2015-09-01-0_RS/data002/data002'] ; % low contrast increments 
DataBlock(1).MbPath{3} = [RootAnalysisDir,'2015-09-01-0_RS/data003/data003'] ; % high contrast decrements 
DataBlock(1).MbPath{4} = [RootAnalysisDir,'2015-09-01-0_RS/data004/data004'] ; % low contrast decrements 
DataBlock(1).DgPath{1} = [RootAnalysisDir,'2015-09-01-0_RS/data005/data005'] ; % lots more of this drifting grating data but only one repeat each 
DataBlock(1).SpontPath = [RootAnalysisDir,'2015-09-01-0_RS/data009/data009'] ; % grey screen
DataBlock(1).FfStepPath{1} = [RootAnalysisDir,'2015-09-01-0_RS/data010/data010'] ; % ff pulses switches light level every 3s
DataBlock(1).mapEiFlag = 'true' ;
DataBlock(1).DSmatFilePath = RootAnalysisSaveDir ;

% Megf10 homozygous KO
% includes image of stb
DataBlock(2).MbPath{1} = [RootAnalysisDir,'2015-09-02-0_RS/data001/data001'] ; % high contrast increments 
DataBlock(2).MbPath{2} = [RootAnalysisDir,'2015-09-02-0_RS/data002/data002'] ; % low contrast increments 
DataBlock(2).MbPath{3} = [RootAnalysisDir,'2015-09-02-0_RS/data003/data003'] ; % high contrast decrements 
DataBlock(2).MbPath{4} = [RootAnalysisDir,'2015-09-02-0_RS/data004/data004'] ; % low contrast decrements 
%DataBlock(2).MbPath{5} = [RootAnalysisDir,'2015-09-02-0_RS/data005/data005'] ; % high contrast decrements 
%DataBlock(2).MbPath{6} = [RootAnalysisDir,'2015-09-02-0_RS/data006/data006'] ; % low contrast decrements 
DataBlock(2).DgPath{1} = [RootAnalysisDir,'2015-09-02-0_RS/data007/data007'] ; % 50% contrast sp 60; tp 1,4; dir 12; 1x
DataBlock(2).DgPath{2} = [RootAnalysisDir,'2015-09-02-0_RS/data008/data008'] ; % 50% contrast sp 60; tp 1,4; dir 12; 1x
DataBlock(2).DgPath{3} = [RootAnalysisDir,'2015-09-02-0_RS/data009/data009'] ; % 50% contrast sp 60; tp 1,4; dir 12; 1x
DataBlock(2).DgPath{4} = [RootAnalysisDir,'2015-09-02-0_RS/data010/data010'] ; % 50% contrast sp 60; tp 1,4; dir 12; 1x
DataBlock(2).SpontPath = [RootAnalysisDir,'2015-09-02-0_RS/data011/data011'] ; % grey screen
DataBlock(2).FfStepPath{1} = [RootAnalysisDir,'2015-09-02-0_RS/data012/data012'] ; % ff pulses switches light level every 3s
DataBlock(2).mapEiFlag = 'true' ;
DataBlock(2).DSmatFilePath = RootAnalysisSaveDir ;

% Megf10 het (previous to stim trigger logic change 2015-09-01)
DataBlock(3).BwPath{1} = [RootAnalysisDir,'2015-06-09-0/data000/data000'] ; % ????
DataBlock(3).MbPath{1} = [RootAnalysisDir,'2015-06-09-0/data001/data001'] ; % high contrast increments 
DataBlock(3).MbPath{2} = [RootAnalysisDir,'2015-06-09-0/data002/data002'] ; % low contrast increments 
DataBlock(3).MbPath{3} = [RootAnalysisDir,'2015-06-09-0/data003/data003'] ; % high contrast decrements 
DataBlock(3).MbPath{4} = [RootAnalysisDir,'2015-06-09-0/data004/data004'] ; % low contrast decrements 
%DataBlock(3).MbPath{5} = [RootAnalysisDir,'2015-06-09-0/data005/data005'] ; % high contrast decrements 
%DataBlock(3).MbPath{6} = [RootAnalysisDir,'2015-06-09-0/data006/data006'] ; % low contrast decrements 
DataBlock(3).DgPath{1} = [RootAnalysisDir,'2015-06-09-0/data007/data007'] ; % 50% contrast sp 60; tp 1,4; dir 12; 1x
DataBlock(3).DgPath{2} = [RootAnalysisDir,'2015-06-09-0/data008/data008'] ; % 50% contrast sp 60; tp 1,4; dir 12; 1x
DataBlock(3).DgPath{3} = [RootAnalysisDir,'2015-06-09-0/data009/data009'] ; % 50% contrast sp 60; tp 1,4; dir 12; 1x
DataBlock(3).DgPath{4} = [RootAnalysisDir,'2015-06-09-0/data010/data010'] ; % 50% contrast sp 60; tp 1,4; dir 12; 1x
DataBlock(3).SpontPath = [RootAnalysisDir,'2015-06-09-0/data011/data011'] ; % grey screen
DataBlock(3).FfStepPath{1} = [RootAnalysisDir,'2015-06-09-0/data012/data012'] ; % ff pulses switches light level every 3s
DataBlock(3).mapEiFlag = 'true' ;
DataBlock(3).DSmatFilePath = RootAnalysisSaveDir ;

% Megf10 het 
DataBlock(4).BwPath{1} = [RootAnalysisDir,'2015-10-09-0_RS/data000/data000'] ; % (10-1-.48 60x60)
DataBlock(4).BwMoviePath{1} = [RootMovieDir,'BW-10-1-0.48-11111-60x60-60.35.xml'] ;
DataBlock(4).MbPath{1} = [RootAnalysisDir,'2015-10-09-0_RS/data001/data001'] ; % high contrast increments 
DataBlock(4).MbPath{2} = [RootAnalysisDir,'2015-10-09-0_RS/data002/data002'] ; % low contrast increments 
DataBlock(4).MbPath{3} = [RootAnalysisDir,'2015-10-09-0_RS/data003/data003'] ; % high contrast decrements 
DataBlock(4).MbPath{4} = [RootAnalysisDir,'2015-10-09-0_RS/data004/data004'] ; % low contrast decrements 
DataBlock(4).DgPath{1} = [RootAnalysisDir,'2015-10-09-0_RS/data005/data005'] ; % 
DataBlock(4).DgPath{2} = [RootAnalysisDir,'2015-10-09-0_RS/data006/data006'] ; % 
DataBlock(4).FfStepPath{1} = [RootAnalysisDir,'2015-06-09-0_RS/data007/data007'] ; % ff pulses switches light level every 3s
DataBlock(4).mapEiFlag = 'true' ;
DataBlock(4).DSmatFilePath = RootAnalysisSaveDir ;

% Wt (piece was small on large array)
DataBlock(5).MbPath{1} = [RootAnalysisDir,'2015-10-22-0_RS/data000/data000'] ; % 100% contrast increments 
DataBlock(5).MbPath{2} = [RootAnalysisDir,'2015-10-22-0_RS/data001/data001'] ; % ~25% contrast increments 
DataBlock(5).MbPath{3} = [RootAnalysisDir,'2015-10-22-0_RS/data002/data002'] ; % 100% contrast decrements
DataBlock(5).MbPath{4} = [RootAnalysisDir,'2015-10-22-0_RS/data003/data003'] ; % ~25% contrast decrements 
DataBlock(5).FfStepPath{1} = [RootAnalysisDir,'2015-10-22-0_RS/data004/data004'] ; % ff pulses switches light level every 3s (less than 25 repeats)
DataBlock(5).mapEiFlag = 'true' ;
DataBlock(5).DSmatFilePath = RootAnalysisSaveDir ;

% TEMP TO TROUBLE SHOOT FUNCTION FROM CNG DATA SET!!
DataBlock(9).BwPath{1} = [RootAnalysisDir,'2015-07-24-0/data005/data005'] ; % 0 NDF (BW 10-1-.48)
DataBlock(9).BwMoviePath{1} = [RootMovieDir,'BW-10-1-0.48-11111-60x60-60.35.xml'] ;

% Params defaults
Params.empty = [] ;
ForIgor = [] ;

% run functions

if perform.MbAnalyzer ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = MbAnalyzer(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.DsCellBwAnalyzer ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = DsCellBwAnalyzer(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end


if perform.MBAnalyzerExtra ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        ForIgor = MBAnalyzerExtra(DataBlock, DB, ForIgor) ;
        %ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

%% export to hdf5 file

% check that field names are not too long for Igor (must be less than 32 characters)
ForIgorFields = fieldnames(ForIgor) ;
for a=1:length(ForIgorFields) ; % for every field
    fieldLength(a) = length(ForIgorFields{a}) ;
end

% repository status NEEDS TO BE FIXED!
temp = getGitInfo ; 
RepVer = temp.hash ; % repository version 
[tempS,tempR] = system('git status') ;
if length(tempR)==62 ;
    ForIgor.RepStat = ['GitHub up to date ',RepVer] ; % if no unsynced files
else
    ForIgor.RepStat = ['GitHub not up to date ',RepVer] ;
end

% if max(fieldLength)<32 ;
%     cd Z:/cafaro' Documents'/Analysis/ForIgorHdfs/   
%     exportStructToHDF5(ForIgor,'PNadaptationAnalysisI.h5','/')
% else disp('name too long for igor')
% end


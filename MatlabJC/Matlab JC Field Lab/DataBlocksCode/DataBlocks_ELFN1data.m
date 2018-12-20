% script to analyze CNG-/- +/- rescue, wt

% 3/17/2015 JC

% functions to perform
perform.BwAnalysis = false ;
perform.BwAnalysisV2 = false ;
perform.BwAnalyzerNoMaping = false ;

perform.DfAnalyzer = false ;
perform.DfAnalyzerV2 = false ;
perform.DfAnalyzerV3 = false ;
perform.DfAnalyzerV4 = false ;
perform.DfAnalyzerV5 = false ;
perform.DimFlashAnalysisWithCellTypes = false ;
perform.DfAnalyzerV3WithBWmapping = false ;

perform.SpontRateAnalysis = false ;
perform.SpontRateAnalysisNoMapping = false ;

perform.DsCellFinder = false ;

perform.KoReverseGratingsAnalyzerConcatPlusMappingMultiSetV1 = true ;

% perform functions on this data block set
DBset_id = 1 ;

% data block sets
DBset{1} = [3] ; % 
DBset{2} = [1,2,3] ; % ELFN1 dim flash

% Root dircetory for analysis files
%RootAnalysisDir = '/Volumes/Berlioz/Analysis/' ; % Berlioz
RootAnalysisDir = '/Volumes/lab/Experiments/Array/Analysis/' ; % Brahms server
RootMovieDir = '/Volumes/lab/acquisition/movie-xml/' ; % Brahms server xml movie root

% ELFN1 from Kirill M Lab (RRL mark)
% sevaral noise events in the dim flash and spont data
% generally good recording
DataBlock(1).SpontPath{1} = [RootAnalysisDir,'2017-05-30-0/data000/data000'] ; % in darkness
DataBlock(1).DfPath = [RootAnalysisDir,'2017-05-30-0/data001/data001'] ; % 3 NDF
DataBlock(1).DfParams.NDF =   [5,5,5,5,5,4,4,4,4,4,4,3,3,3,3,3,3,2,2,2,2,2,2,1,1,1,1,1,1] ; % on filter turret
DataBlock(1).DfParams.Ftime = [8,4,2,8,4,4,2,8,4,2,8,8,4,2,4,8,2,2,4,8,2,4,8,8,2,4,2,8,4] ; % ms
DataBlock(1).DfParams.interFlashInt = 3 ; % sec
DataBlock(1).BwPath{1} = [RootAnalysisDir,'2017-05-30-0/data005/data005'] ; % 2 NDF (BW 30-4-.48)
DataBlock(1).BwPath{2} = [RootAnalysisDir,'2017-05-30-0/data006/data006'] ; % 0 NDF (BW 10-1-.48)
DataBlock(1).DfBwBlock = 2 ; % BW block to use for classificiation 
DataBlock(1).ExampleCelli = 1 ; % 
DataBlock(1).DsPath = [RootAnalysisDir,'2017-05-30-0/data007/data007'] ; % 0 NDF (2 TP, 1 SP)
DataBlock(1).FfStepPath{1} = [RootAnalysisDir,'2017-05-30-0/data003/data003'] ; % 3 NDF gwgb

% ELFN1 from Kirill M Lab (cage KM 647 L ear mark)
% not a great recording
DataBlock(2).DfPath = [RootAnalysisDir,'2017-06-28-0/data000/data000'] ; % 3 NDF
DataBlock(2).DfParams.NDF =   [5,4,3,2,2,2,2,2,2,1,1,1,1,1,1,0,0,0] ; % on filter turret
DataBlock(2).DfParams.Ftime = [8,8,8,8,2,4,2,8,4,2,4,2,8,4,8,2,4,2] ; % ms
DataBlock(2).DfParams.interFlashInt = 3 ; % sec
DataBlock(2).BwPath{1} = [RootAnalysisDir,'2017-06-28-0/data001/data001'] ; % 0 NDF (BW 15-2)
DataBlock(2).DfBwBlock = 1 ; % BW block to use for classificiation 
DataBlock(2).ExampleCelli = 1 ; % 
DataBlock(2).DsPath = [RootAnalysisDir,'2017-06-28-0/data002/data002'] ; % 0 NDF (2 TP, 1 SP)

% ELFN1 from Kirill M Lab (cage 647 RR mark)
% generally good recording
DataBlock(3).SpontPath{1} = [RootAnalysisDir,'2017-07-12-0/data000-1/data000-1'] ; % in darkness
DataBlock(3).DfPath = [RootAnalysisDir,'2017-07-12-0/data002/data002'] ; % 3 NDF
DataBlock(3).DfParams.NDF =   [5,4,3,3,3,3,2,2,2,2,2,2,1,1,1,1,1,1,0] ; % on filter turret
DataBlock(3).DfParams.Ftime = [8,8,8,4,8,4,2,4,2,8,4,8,2,4,8,2,4,8,2] ; % ms
DataBlock(3).DfParams.interFlashInt = 3 ; % sec
DataBlock(3).BwPath{1} = [RootAnalysisDir,'2017-07-12-0/data003/data003'] ; % 2 NDF (BW 15-2-.48)
DataBlock(3).BwPath{2} = [RootAnalysisDir,'2017-07-12-0/data004/data004'] ; % 2 NDF (BW 30-4-.48)
DataBlock(3).BwPath{3} = [RootAnalysisDir,'2017-07-12-0/data005/data005'] ; % 0 NDF (BW 10-1-.48) - no ON cells
DataBlock(3).DfBwBlock = 2 ; % BW block to use for classificiation 
DataBlock(3).ExampleCelli = 1 ; % 
DataBlock(3).DsPath = [RootAnalysisDir,'2017-07-12-0/data006/data006'] ; % 0 NDF (2 TP, 1 SP)
DataBlock(3).RgPath{1} = [RootAnalysisDir,'2017-07-12-0/data007/data007'] ; % 0 NDF reverse grating 50% contrast 3 phases, 1 TP, 6 SP
DataBlock(3).RgPath{1} = [RootAnalysisDir,'2017-07-12-0/data008/data008'] ; % 0 NDF reverse grating 100% contrast 3 phases, 1 TP, 6 SP
DataBlock(3).RepStim{1} = [RootAnalysisDir,'2017-07-12-0/data009/data009'] ; % repeated cat movie

% Params defaults
Params.empty = [] ;
ForIgor = [] ;

% run functions

if perform.BwAnalysis ;
    for a=1:length(DBset(DBset_id)) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = BwAnalyzer(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.BwAnalysisV2 ;
    for a=1:length(DBset(DBset_id)) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = BwAnalyzerV2(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.BwAnalyzerNoMaping ;
    for a=1:length(DBset(DBset_id)) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = BwAnalyzerNoMaping(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.DfAnalyzer ;
    for a=1:length(DBset(DBset_id)) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = DfAnalyzer(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end    

if perform.DfAnalyzerV2 ;
    for a=1:length(DBset(DBset_id)) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = DfAnalyzerV2(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end    

if perform.DfAnalyzerV3 ;
    for a=1:length(DBset(DBset_id)) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = DfAnalyzerV3(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end  

if perform.DfAnalyzerV4 ;
    for a=1:length(DBset(DBset_id)) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = DfAnalyzerV4(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end  

if perform.DfAnalyzerV5 ;
    for a=1:length(DBset(DBset_id)) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = DfAnalyzerV5(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end 

if perform.DfAnalyzerV3WithBWmapping ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = DfAnalyzerV3WithBWmapping(DataBlock, DB, Params) ;
        ForMatlab(DB) = Temp.ForMatlab(DB) ; % put into structure for matlab
        Temp = rmfield(Temp,'ForMatlab') ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end  

if perform.DimFlashAnalysisWithCellTypes ;
    for a=1:length(DBset(DBset_id)) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = DimFlashAnalysisWithCellTypes(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end    

if perform.SpontRateAnalysis ;
    for a=1:length(DBset(DBset_id)) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = SpontRateAnalysis(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end  

if perform.SpontRateAnalysisNoMapping ;
    for a=1:length(DBset(DBset_id)) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = SpontRateAnalysisNoMapping(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end 


if perform.DsCellFinder ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = DsCellFinder(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end


if perform.KoReverseGratingsAnalyzerConcatPlusMappingMultiSetV1  ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = KoReverseGratingsAnalyzerConcatPlusMappingMultiSetV1 (DataBlock, DB, Params) ;
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


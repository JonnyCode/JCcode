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

% perform functions on this data block set
DBset_id = 12 ;

% data block sets
DBset{1} = [32] ; % 
DBset{2} = [1,4,5,13,14,19,20,21,30,31] ; % CNG KO 
DBset{3} = [2,3,6,15,16] ; % CNG KO rescue
DBset{4} = [7,8,9] ; % wt mouse
DBset{5} = [18] ; % Cre-Er mouse
DBset{6} = [12,25,27] ; % CNG deltaCam
DBset{7} = [13,19,20,21,30,31] ; % dim flash (modern protocol) CNG KO at 3m
DBset{8} = [3,23,24] ; % dim flash (modern protocol) CNG rescue at 1m
DBset{9} = [22] ; % dim flash (modern protocol) CNGB - Cre-negative rescue at 1m
DBset{10} = [26,28,29] ; % dim flash (modern protocol) CNG KO at 1m
DBset{11} = [30,31,3,23,24,12,25,27,19,20,21,26,28,29] ; % CNG for paper
DBset{12} = [32,33,34] ; % ELFN1

% Root dircetory for analysis files
%RootAnalysisDir = '/Volumes/Berlioz/Analysis/' ; % Berlioz
RootAnalysisDir = '/Volumes/lab/Experiments/Array/Analysis/' ; % Brahms server
RootMovieDir = '/Volumes/lab/acquisition/movie-xml/' ; % Brahms server xml movie root

% CNG KO from Jeannie Chen's lab
% XyY analyzed dim flash from this mouse
% 80/20 NDF 2 at beginging then switched to 50/50 0 NDF before data002 
DataBlock(1).BwPath{1} = [RootAnalysisDir,'2014-01-14-0/data008/data008'] ; % 0 NDF
DataBlock(1).DfPath = [RootAnalysisDir, '2014-01-14-0/data000-1/data000-1'] ; % 2 NDF 
DataBlock(1).DfParams.NDF = [4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5] ;
DataBlock(1).DfParams.Ftime = [1,1,1,1,10,10,10,10,10,10,20,50,50,100,10,10,10,1,1,10,10] ; % ms
DataBlock(1).DfParams.Setting = [1,.6,.6,1,1,.8,.6,.7,1,3,3,3,3,3,1,.6,.8,1,.6,.7,1] ;
DataBlock(1).DfParams.interFlashInt = 3 ; % sec
DataBlock(1).FfStepPath{1} = [RootAnalysisDir, '2014-01-14-0/data002/data002'] ; % 5 NDF wgbg 
DataBlock(1).FfStepPath{2} = [RootAnalysisDir, '2014-01-14-0/data003/data003'] ; % 4 NDF wgbg 
DataBlock(1).FfStepPath{3} = [RootAnalysisDir, '2014-01-14-0/data004/data004'] ; % 3 NDF wgbg 
DataBlock(1).FfStepPath{4} = [RootAnalysisDir, '2014-01-14-0/data005/data005'] ; % 2 NDF wgbg 
DataBlock(1).FfStepPath{5} = [RootAnalysisDir, '2014-01-14-0/data006/data006'] ; % 1 NDF wgbg 
DataBlock(1).FfStepPath{6} = [RootAnalysisDir, '2014-01-14-0/data007/data007'] ; % 0 NDF wgbg 
DataBlock(1).DsPath = [RootAnalysisDir,'2014-01-14-0/data009/data009'] ; % 0 NDF (2 TP, 1 SP, 8 directions)
DataBlock(1).mapEiFlag = 'true' ;

% CNG KO rescue from Jeannie Chen's lab
% XyY analyzed dim flash from this mouse
% 80/20 NDF 2 at beginging (then switched to 50/50 0 NDF before data002 ????????) 
DataBlock(2).BwPath{1} = [RootAnalysisDir,'2014-01-16-0/data008/data008'] ; % 0 NDF
DataBlock(2).DfPath = [RootAnalysisDir, '2014-01-16-0/data000-1/data000-1'] ; % 2 NDF 
DataBlock(2).DfParams.NDF = [5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6] ;
DataBlock(2).DfParams.Ftime = [1,1,1,1,10,10,10,10,1,1,1,10,10,10,10,1,1,1,10,10,20,50,10] ; % ms
DataBlock(2).DfParams.Setting = [1,.6,1,.6,.6,1,.8,.6,3,2,1,.6,.7,1.2,.8,1.1,.8,.6,1,3,3,3,.6] ;
DataBlock(2).DfParams.interFlashInt = 3 ; % sec
DataBlock(2).FfStepPath{1} = [RootAnalysisDir, '2014-01-16-0/data002/data002'] ; % 5 NDF wgbg 
DataBlock(2).FfStepPath{2} = [RootAnalysisDir, '2014-01-16-0/data003/data003'] ; % 4 NDF wgbg 
DataBlock(2).FfStepPath{3} = [RootAnalysisDir, '2014-01-16-0/data004/data004'] ; % 3 NDF wgbg 
DataBlock(2).FfStepPath{4} = [RootAnalysisDir, '2014-01-16-0/data005/data005'] ; % 2 NDF wgbg 
DataBlock(2).FfStepPath{5} = [RootAnalysisDir, '2014-01-16-0/data006/data006'] ; % 1 NDF wgbg 
DataBlock(2).FfStepPath{6} = [RootAnalysisDir, '2014-01-16-0/data007/data007'] ; % 0 NDF wgbg 
DataBlock(2).DsPath = [RootAnalysisDir,'2014-01-16-0/data009/data009'] ; % 0 NDF (2 TP, 1 SP, 8 directions)
DataBlock(2).mapEiFlag = 'true' ;

% CNG KO rescue (TM at 1m record at ~5m) from Jeannie Chen's lab
% "mapped array and took image of the retina"
% 80/20 NDF 3 entire time but then switch to 50/50 NDF 0 before data002
DataBlock(3).BwPath{1} = [RootAnalysisDir,'2014-10-24-0_JC/data005/data005'] ; % 0 NDF (BW 10-1-.48 3600s)
DataBlock(3).BwPath{2} = [RootAnalysisDir,'2014-10-24-0_JC/data002/data002'] ; % 3 NDF (BW 10-4-.48 3600s)
DataBlock(3).BwMoviePath{1} = [RootMovieDir,'BW-10-1-0.48-11111-60x60-60.35.xml'] ;
DataBlock(3).BwMoviePath{2} = [RootMovieDir,'BW-10-4-0.48-11111-60x60-60.35.xml'] ;
DataBlock(3).DfPath = [RootAnalysisDir,'2014-10-24-0_RS/data001-RS/data001-RS'] ; % 3 NDF
DataBlock(3).DfParams.NDF =   [5,5,4,5,4,5,5,4,3,4,4,3,3,2,4,2,1,0,4] ; % on filter turret
DataBlock(3).DfParams.Ftime = [2,4,2,8,4,4,2,2,2,4,8,4,8,2,2,8,2,2,2] ; % ms
DataBlock(3).DfParams.interFlashInt = 3 ; % sec
DataBlock(3).SpontPath{1} = [RootAnalysisDir,'2014-10-24-0_JC/data000/data000'] ; % in darkness
DataBlock(3).FfStepPath{1} = [RootAnalysisDir,'2014-10-24-0_JC/data003/data003'] ; % 3 NDF gwgb
DataBlock(3).FfStepPath{2} = [RootAnalysisDir,'2014-10-24-0_JC/data009/data009'] ; % 0 NDF gwgb
DataBlock(3).BwRepeatPath{1} = [RootAnalysisDir,'2014-10-24-0_JC/data007/data007'] ; % 0 NDF (BW 10-1-.48 120 repeats 10 s epochs) 
DataBlock(3).DsPath = [RootAnalysisDir,'2014-10-24-0_JC/data008/data008'] ; % 0 NDF (2 TP, 2 SP, 18 directions)
DataBlock(3).mapEiFlag = 'true' ;
DataBlock(3).DfBwBlock = 1 ; % BW block to use for classificiation 
DataBlock(3).ExampleCelli = 1 ; % random cell because there are no matches

% CNG KO from Jeannie Chen's lab (recorded around 1m)
% Notes issue with data001 (16 df blocks noted - assumed NDF1, 2ms was a
% duplicate and not delivered)
% 80/20 NDF 3 start but then switch to 50/50 NDF 0 before data002
DataBlock(4).BwPath{1} = [RootAnalysisDir,'2014-10-29-0_JC/data007/data007'] ; % 0 NDF (BW  10-2-.48) 3600 s)
DataBlock(4).BwPath{2} = [RootAnalysisDir,'2014-10-29-0_JC/data005/data005'] ; % 2 NDF (BW 20-6-.48 1800 s)
DataBlock(4).BwPath{3} = [RootAnalysisDir,'2014-10-29-0_JC/data003/data003'] ; % 3 NDF (BW 30-6-.48 600s)
DataBlock(4).BwPath{4} = [RootAnalysisDir,'2014-10-29-0_JC/data002/data002'] ; % 3 NDF (BW 10-4-.48 600 s)
DataBlock(4).BwMoviePath{1} = [RootMovieDir,'BW-10-2-0.48-11111-60x60-60.35.xml'] ;
DataBlock(4).BwMoviePath{2} = [RootMovieDir,'BW-20-6-0.48-11111-30x30-60.35.xml'] ;
DataBlock(4).BwMoviePath{3} = [RootMovieDir,'BW-30-6-0.48-11111-20x20-60.35.xml'] ;
DataBlock(4).BwMoviePath{4} = [RootMovieDir,'BW-10-4-0.48-11111-60x60-60.35.xml'] ;
DataBlock(4).DfPath = [RootAnalysisDir,'2014-10-29-0_RS/data001-RS/data001-RS'] ; % 3 NDF
DataBlock(4).DfParams.NDF =   [5,4,4,3,2,3,3,3,3,2,2,1,2,1,0] ; % on filter turret
DataBlock(4).DfParams.Ftime = [2,2,8,2,2,4,2,8,4,4,2,2,8,4,1] ; % ms
DataBlock(4).DfParams.interFlashInt = 3 ; % sec
DataBlock(4).SpontPath{1} = [RootAnalysisDir,'2014-10-29-0_JC/data000/data000'] ; % in darkness
DataBlock(4).FfStepPath{1} = [RootAnalysisDir,'2014-10-29-0_JC/data004/data004'] ; % 3 NDF gwgb
DataBlock(4).FfStepPath{2} = [RootAnalysisDir,'2014-10-29-0_JC/data006/data006'] ; % 2 NDF gwgb
DataBlock(4).FfStepPath{3} = [RootAnalysisDir,'2014-10-29-0_JC/data008/data008'] ; % 0 NDF gwgb
DataBlock(4).BwRepeatPath{1} = [RootAnalysisDir,'2014-10-29-0_JC/data010/data010'] ; % 0 NDF (BW 10-1-.48 120 repeats 10 s epochs) 
DataBlock(4).DsPath = [RootAnalysisDir,'2014-10-29-0_JC/data009/data009'] ; % 0 NDF (18 directions)
DataBlock(4).mapEiFlag = 'true' ;

% CNG KO from Jeannie Chen's lab
% 80/20 NDF 3 start (but then switch to 50/50 NDF 0 before data002 ????)
DataBlock(5).BwPath{1} = [RootAnalysisDir,'2013-12-17-0/data006/data006'] ; % 0 NDF (BW 15-1-.48-11111 40x40)
DataBlock(5).BwPath{2} = [RootAnalysisDir,'2013-12-17-0/data010/data010'] ; % 0 NDF (BW 1-8-.48-11111 600x600)
DataBlock(5).BwPath{3} = [RootAnalysisDir,'2013-12-17-0/data011/data011'] ; % 0 NDF (BW 15-1-.48-11111 40x40)
DataBlock(5).DfPath = [RootAnalysisDir, '2013-12-17-0/data000-1/data000-1'] ; % 2 NDF 
DataBlock(5).DfParams.NDF = [4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4] ;
DataBlock(5).DfParams.Ftime = [1,1,1,10,10,10,20,20,20,20,50,50,50,50,10,10,10,20,1,10,1,50] ; % ms
DataBlock(5).DfParams.Setting = [3,1,5,1,3,5,1,3,2,5,1,2,3,5,3,1,5,2,5,.7,.7,3] ;
DataBlock(5).DfParams.interFlashInt = 3 ; % sec
DataBlock(5).FfStepPath{1} = [RootAnalysisDir, '2013-12-17-0/data002/data002'] ; % 2 NDF wgbg 
DataBlock(5).FfStepPath{2} = [RootAnalysisDir, '2013-12-17-0/data003/data003'] ; % 3 NDF wgbg 
DataBlock(5).FfStepPath{3} = [RootAnalysisDir, '2013-12-17-0/data004/data004'] ; % 4 NDF wgbg 
DataBlock(5).FfStepPath{4} = [RootAnalysisDir, '2013-12-17-0/data007/data007'] ; % 0 NDF wgbg 
DataBlock(5).DsPath = [RootAnalysisDir,'22013-12-17-0/data009/data009'] ; % 0 NDF (2 TP, 1 SP, 8 directions)
DataBlock(5).mapEiFlag = 'true' ;

% CNG KO rescue from Jeannie Chen's lab
% 80/20 NDF 3 start but then switch to 50/50 NDF 0 before data00
DataBlock(6).BwPath{1} = [RootAnalysisDir,'2013-12-20-0/data004/data004'] ; % 0 NDF (10-1-.48 60x60)
DataBlock(6).DfPath = [RootAnalysisDir, '2013-12-20-0/data000/data000'] ; % 2 NDF 
DataBlock(6).DfParams.NDF = [5,5,5,5,5,5,5,5,5,5,5] ;
DataBlock(6).DfParams.Ftime = [1,10,10,20,1,1,1,10,10,10,10] ; % ms
DataBlock(6).DfParams.Setting = [1,1,2,2,2,4,.6,1,.8,.7,.6] ;
DataBlock(6).DfParams.interFlashInt = 3 ; % sec
DataBlock(6).FfStepPath{1} = [RootAnalysisDir, '2013-12-20-0/data002/data002'] ; % 5 NDF wgbg 
DataBlock(6).FfStepPath{2} = [RootAnalysisDir, '2013-12-20-0/data003/data003'] ; % 4 NDF wgbg 
DataBlock(6).FfStepPath{3} = [RootAnalysisDir, '2013-12-20-0/data005/data005'] ; % 0 NDF wgbg 
DataBlock(6).DsPath = [RootAnalysisDir,'22013-12-20-0/data006/data006'] ; % 0 NDF (2 TP, 1 SP, 8 directions)
DataBlock(6).mapEiFlag = 'true' ;

% wt mouse c57bl/6 (DOB - 4/14/2015)
% 80/20 NDF 3 but than switched to 50/50 NDF 0 before data004
% spikes on left of array became apparent around data009 - not sure why
DataBlock(7).BwPath{1} = [RootAnalysisDir,'2015-06-17-0/data009/data009'] ; % 0 NDF (BW 10-1-.48 3600s)
DataBlock(7).BwPath{2} = [RootAnalysisDir,'2015-06-17-0/data004/data004'] ; % 2 NDF (BW 10-2-.48 3600s)
DataBlock(7).BwPath{3} = [RootAnalysisDir,'2015-06-17-0/data006/data006'] ; % 2 NDF (BW 30-4-.48 3600s)
DataBlock(7).BwMoviePath{1} = [RootMovieDir,] ;
DataBlock(7).BwMoviePath ;
DataBlock(7).DfPath = [RootAnalysisDir,'2015-06-17-0/data002-3/data002-3'] ; % 3 NDF
DataBlock(7).DfParams.NDF =   [4,4,4,4,5,5,5,4,4,4,4,4,4,5,3,3,3,3,3,2,2,2,1,1] ; % on filter turret (guess to correct error in notes)
DataBlock(7).DfParams.Ftime = [2,2,8,4,4,2,8,8,2,2,4,2,8,8,2,4,8,8,2,2,8,4,4,8] ; % ms
% DataBlock(7).DfParams.NDF =   [4,4,4,4,4,5,5,5,4,4,4,4,4,4,5,3,3,3,3,3,2,2,2,1,1] ; % on filter turret
% DataBlock(7).DfParams.Ftime = [2,2,2,8,4,4,2,8,8,2,2,4,2,8,8,2,4,8,8,2,2,8,4,4,8] ; % ms
DataBlock(7).DfParams.interFlashInt = [3,5] ; % sec
DataBlock(7).SpontPath{1} = [RootAnalysisDir,'2015-06-17-0/data001/data001'] ; % in darkness
DataBlock(7).SpontPath{2} = [RootAnalysisDir,'2015-06-17-0/data010/data010'] ; % NDF 0 OLED gray screen
DataBlock(7).FfStepPath{1} = [RootAnalysisDir,'2015-06-17-0/data008/data008'] ; % 2 NDF gwgb
DataBlock(7).FfStepPath{2} = [RootAnalysisDir,'2015-06-17-0/data013/data013'] ; % 0 NDF gwgb
DataBlock(7).BwRepeatPath{1} = [RootAnalysisDir,'2015-06-17-0/data005/data005'] ; % 2 NDF (BW 10-2-.48) 
DataBlock(7).BwRepeatPath{2} = [RootAnalysisDir,'2015-06-17-0/data007/data007'] ; % 2 NDF (BW 30-4-.48) 
DataBlock(7).BwRepeatPath{3} = [RootAnalysisDir,'2015-06-17-0/data011/data011'] ; % 0 NDF (BW 10-1-.48) 
DataBlock(7).DsPath = [RootAnalysisDir,'2015-06-17-0/data012/data012'] ; % 0 NDF (2 TP, 1 SP, 12 directions)
DataBlock(7).mapEiFlag = 'true' ;

% wt mouse c57bl/6 (DOB - 6/12/2015)
% 80/20 NDF 3 
DataBlock(8).DfPath = [RootAnalysisDir,'2015-07-09-0/data002/data002'] ; % 3 NDF
DataBlock(8).DfParams.NDF =   [2,2,5,5,5,5,5,5,4,2,4,4,4,4,4,5,2,3,3,3] ; % NDF (guess to correct error in notes)
DataBlock(8).DfParams.Ftime = [8,2,2,8,4,8,2,4,8,2,2,4,2,8,4,4,2,2,4,8] ; % ms
DataBlock(8).DfParams.interFlashInt = [3] ; % sec
DataBlock(8).mapEiFlag = 'true' ;

% wt mouse c57bl/6 (DOB - 4/14/2015)
% 80/20 NDF 3 but than switched to 50/50 NDF 0 before data003
% no RPE on this retina
DataBlock(9).BwPath{1} = [RootAnalysisDir,'2015-07-24-0/data005/data005'] ; % 0 NDF (BW 10-1-.48)
DataBlock(9).BwPath{2} = [RootAnalysisDir,'2015-07-24-0/data006/data006'] ; % 0 NDF (BW 30-4-.48)
DataBlock(9).BwPath{3} = [RootAnalysisDir,'2015-07-24-0/data003/data003'] ; % 5? NDF (BW 10-2-.48 3600s)
DataBlock(9).BwPath{4} = [RootAnalysisDir,'2015-07-24-0/data004/data004'] ; % 5? NDF (BW 30-4-.48 3600s)
DataBlock(9).BwMoviePath{1} = [RootMovieDir,'BW-10-1-0.48-11111-60x60-60.35.xml'] ;
%DataBlock(9).BwMoviePath{2} ;
%DataBlock(9).DfPath = [RootAnalysisDir,'2015-07-24-0/data001/data001'] ; % 3 NDF (this data has some issues - easiest to avoid
DataBlock(9).DfPath = [RootAnalysisDir,'2015-07-24-0/data002/data002'] ; % 3 NDF
DataBlock(9).DfParams.NDF =   [5,5,5,4,4,4,2,2,2,3,0,1] ; % on filter turret 
DataBlock(9).DfParams.Ftime = [4,8,2,2,4,8,2,4,8,2,2,2] ; % ms
DataBlock(9).DfParams.interFlashInt = [3] ; % sec
DataBlock(9).SpontPath{1} = [RootAnalysisDir,'2015-07-24-0/data000/data000'] ; % in darkness
DataBlock(9).FfStepPath{1} = [RootAnalysisDir,'2015-07-24-0/data008/data008'] ; % 0 NDF gwgb
DataBlock(9).DsPath = [RootAnalysisDir,'2015-07-24-0/data009/data009'] ; % 0 NDF (2 TP, 1 SP, 12 directions)
DataBlock(9).mapEiFlag = 'true' ;

% wt CBA mouse (DOB 2009-06-29)
% 80/20 NDF 3 but than switched to 50/50 NDF 0 before data002
% no RPE on this retina
% 10X was used duing dim flashes! 4x for all else (THUS THE APPARENT SENSITIVITY WILL BE LOWER)
DataBlock(10).BwPath{1} = [RootAnalysisDir,'2015-08-25-0/data004/data004'] ; % 0 NDF (BW 10-1-.48)
DataBlock(10).BwPath{2} = [RootAnalysisDir,'2015-08-25-0/data005/data005'] ; % 0 NDF (BW 30-4-.48)
DataBlock(10).BwPath{3} = [RootAnalysisDir,'2015-08-25-0/data002/data002'] ; % 2 NDF (BW 10-2-.48 3600s)
DataBlock(10).BwPath{4} = [RootAnalysisDir,'2015-08-25-0/data003/data003'] ; % 2 NDF (BW 30-4-.48 3600s)
DataBlock(10).BwMoviePath{1} = [RootMovieDir,] ;
DataBlock(10).BwMoviePath ;
DataBlock(10).DfPath = [RootAnalysisDir,'2015-08-25-0/data001/data001'] ; % 3 NDF
DataBlock(10).DfParams.NDF =   [5,5,5,5,5,5,4,4,4,4,4,4,3,3,3,3,3,3,2,2,2,1,0] ; % on filter turret 
DataBlock(10).DfParams.Ftime = [8,2,4,2,4,8,2,4,8,4,2,8,2,4,8,2,4,8,2,4,8,2,2] ; % ms
DataBlock(10).DfParams.interFlashInt = [3] ; % sec
DataBlock(10).SpontPath{1} = [RootAnalysisDir,'2015-08-25-0/data000/data000'] ; % in darkness
DataBlock(10).FfStepPath{1} = [RootAnalysisDir,'2015-08-25-0/data006/data006'] ; % 0 NDF gwgb
DataBlock(10).DsPath = [RootAnalysisDir,'2015-08-25-0/data007/data007'] ; % 0 NDF (2 TP, 2 SP, 12 directions)
DataBlock(10).mapEiFlag = 'true' ;

% CNGdeltaCaM (2014-10-09 experiment)
% has some dim flash data with missing triggers
DataBlock(11).BwPath{1} = [RootAnalysisDir,'2014-10-09-0/data004/data004'] ; % 0 NDF (BW 15-1-.48)
DataBlock(11).BwPath{2} = [RootAnalysisDir,'2014-10-09-0/data001/data001'] ; % 3 NDF (BW 15-4-.48)
DataBlock(11).BwMoviePath{1} = [RootMovieDir,] ;
DataBlock(11).BwMoviePath ;
DataBlock(11).FfStepPath{1} = [RootAnalysisDir,'2014-10-09-0/data002/data002'] ; % 0 NDF gwgb
DataBlock(11).DsPath = [RootAnalysisDir,'2014-08-25-0/data006/data006'] ; % 0 NDF (2 TP, 2 SP, 12 directions)
DataBlock(11).mapEiFlag = 'true' ;

% CNGdeltaCaM (2014-09-23 experiment)
% some ffpulse +APB +DNQX +UBP
DataBlock(12).BwPath{1} = [RootAnalysisDir,'2014-09-23-0/data004/data004'] ; % 0 NDF (BW 10-1-.48)
DataBlock(12).BwPath{2} = [RootAnalysisDir,'2014-09-23-0/data001/data001'] ; % 3 NDF (BW 15-4-.48 3600s)
DataBlock(12).BwMoviePath{1} = [RootMovieDir,] ;
DataBlock(12).BwMoviePath ;
DataBlock(12).DfPath = [RootAnalysisDir,'2014-09-23-0/data000/data000'] ; % 3 NDF
DataBlock(12).DfParams.NDF =   [5,5,5,4,5,5,4,4,5,4,3,4,4,3,3,2,1] ; % on filter turret 
DataBlock(12).DfParams.Ftime = [2,8,4,2,2,4,4,8,8,2,2,4,8,8,4,2,2] ; % ms
DataBlock(12).DfParams.interFlashInt = [3,6] ; % sec
DataBlock(12).FfStepPath{1} = [RootAnalysisDir,'2014-09-23-0/data002/data002'] ; % 3 NDF gwgb
DataBlock(12).FfStepPath{2} = [RootAnalysisDir,'2014-09-23-0/data005/data005'] ; % 0 NDF gwgb
DataBlock(12).DsPath = [RootAnalysisDir,'2014-09-23-0/data006/data006'] ; % 0 NDF (2 TP, 2 SP, 12 directions)
DataBlock(12).mapEiFlag = 'true' ;
DataBlock(12).DfBwBlock = 1 ; % BW block to use for classificiation 
DataBlock(12).ExampleCelli = 1 ; % random cell because there are no matches

% CNG KO (~3months old) from Jeannie Chen's lab
% something wrong trying to get STA datasets?????
DataBlock(13).SpontPath{1} = [RootAnalysisDir,'2015-10-28-0/data000/data000'] ; % in darkness
DataBlock(13).DfPath = [RootAnalysisDir,'2015-10-28-0/data001-2/data001-2'] ; % 3 NDF
DataBlock(13).DfParams.NDF =   [5,5,5,4,4,4,4,4,4,3,3,3,2,2,2,1,1,1,1,0,0,0] ; % on filter turret 
DataBlock(13).DfParams.Ftime = [2,4,8,8,2,4,2,4,8,8,2,4,4,2,8,2,4,4,8,2,4,8] ; % ms
DataBlock(13).DfParams.interFlashInt = [3] ; % sec
DataBlock(13).BwPath{1} = [RootAnalysisDir,'2015-10-28-0/data006/data006'] ; % 0 NDF (BW 10-1-.48)
DataBlock(13).BwPath{2} = [RootAnalysisDir,'2015-10-28-0/data007/data007'] ; % 0 NDF (BW 30-4-.48)
DataBlock(13).BwPath{3} = [RootAnalysisDir,'2015-10-28-0/data004/data004'] ; % 2 NDF (BW 10-2-.48 3600s)
DataBlock(13).BwPath{4} = [RootAnalysisDir,'2015-10-28-0/data005/data005'] ; % 2 NDF (BW 30-4-.48 3600s)
DataBlock(13).BwMoviePath{1} = [RootMovieDir,] ;
DataBlock(13).FfStepPath{1} = [RootAnalysisDir,'2015-10-28-0/data003/data003'] ; % 3 NDF gwgb
DataBlock(13).FfStepPath{1} = [RootAnalysisDir,'2015-10-28-0/data008/data008'] ; % 3 NDF gwgb
DataBlock(13).DsPath = [RootAnalysisDir,'2015-10-28-0/data007/data007'] ; % 0 NDF (2 TP, 2 SP, 12 directions)
DataBlock(13).mapEiFlag = 'true' ;
DataBlock(13).DfBwBlock = 1 ; % BW block to use for classificiation 
DataBlock(13).ExampleCelli = 1 ; % random

% CNG KO (~3months old) from Jeannie Chen's lab
% missing triggers for data000, 001, ... 
DataBlock(14).SpontPath{1} = [RootAnalysisDir,'2015-10-29-0/data000/data000'] ; % in darkness
DataBlock(14).SpontPath{2} = [RootAnalysisDir,'2015-10-29-0/data005/data005'] ; % ndf 0 grey sceen
DataBlock(14).DfPath = [RootAnalysisDir,'2015-10-29-0/data001/data001'] ; % 3 NDF
DataBlock(14).DfParams.NDF =   [5,5,5,4,4,4,4,4,4,3,3,3,3,2,2,2,1,1,1,1] ; % on filter turret 
DataBlock(14).DfParams.Ftime = [2,4,8,2,4,8,2,4,8,2,4,8,2,2,4,8,2,4,8,2] ; % ms
DataBlock(14).DfParams.interFlashInt = [3] ; % sec
DataBlock(14).BwPath{1} = [RootAnalysisDir,'2015-10-29-0/data006/data006'] ; % 0 NDF (BW 10-1-.48) seed 11111
DataBlock(14).BwPath{1} = [RootAnalysisDir,'2015-10-29-0/data008/data008'] ; % 0 NDF (BW 10-1-.48) seed 22222
DataBlock(14).BwPath{2} = [RootAnalysisDir,'2015-10-29-0/data009/data009'] ; % 0 NDF (BW 30-4-.48) COULD BE data011
DataBlock(14).BwPath{3} = [RootAnalysisDir,'2015-10-29-0/data003/data003'] ; % 2 NDF (BW 10-2-.48 3600s)
DataBlock(14).BwPath{4} = [RootAnalysisDir,'2015-10-29-0/data004/data004'] ; % 2 NDF (BW 30-4-.48 3600s)
DataBlock(14).BwMoviePath{1} = [RootMovieDir,] ;
DataBlock(14).FfStepPath{1} = [RootAnalysisDir,'2015-10-29-0/data002/data002'] ; % 3 NDF gwgb
DataBlock(14).FfStepPath{1} = [RootAnalysisDir,'2015-10-29-0/data011/data011'] ; % 0 NDF gwgb COULD BE data009
DataBlock(14).DsPath = [RootAnalysisDir,'2015-10-29-0/data010/data010'] ; % 0 NDF 
DataBlock(14).mapEiFlag = 'true' ;


% CNG KO resecue (TX at 1 m) from Jeannie Chen's lab
% 80/20 NDF 3 but than switched to 50/50 NDF 0 before data002
DataBlock(15).SpontPath{1} = [RootAnalysisDir,'2015-10-30-0/data000/data000'] ; % in darkness
DataBlock(15).DfPath = [RootAnalysisDir,'2015-10-30-0/data001/data001'] ; % 3 NDF
DataBlock(15).DfParams.NDF =   [4,4,4,4,4,3,3,3,2,2,2,1,1,1] ; % on filter turret (UNCLEAR IF THIS IS CORRECT! - assuming the only the last 14 set triggs were recorded)
DataBlock(15).DfParams.Ftime = [4,8,2,4,8,2,4,8,2,4,8,2,4,8] ; % ms
% DataBlock(15).DfParams.NDF =   [5,5,5,5,5,4,4,4,4,4,4,3,3,3,2,2,2,1,1,1] ; % on filter turret
% DataBlock(15).DfParams.Ftime = [2,4,8,4,8,2,4,8,2,4,8,2,4,8,2,4,8,2,4,8] ; % ms
DataBlock(15).DfParams.interFlashInt = 3 ; % sec
DataBlock(15).FfStepPath{1} = [RootAnalysisDir,'2015-10-30-0/data002/data002'] ; % 3 NDF gwgb
DataBlock(15).FfStepPath{2} = [RootAnalysisDir,'2015-10-30-0/data005/data005'] ; % 0 NDF gwgb
DataBlock(15).BwPath{1} = [RootAnalysisDir,'2015-10-30-0/data003/data003'] ; % 2 NDF (BW 10-2-.48)
DataBlock(15).BwPath{2} = [RootAnalysisDir,'2015-10-30-0/data004/data004'] ; % 2 NDF (BW 30-4-.48)
DataBlock(15).BwPath{3} = [RootAnalysisDir,'2015-10-30-0/data006/data006'] ; % 0 NDF (BW 10-1-.48)
DataBlock(15).BwPath{4} = [RootAnalysisDir,'2015-10-30-0/data008/data008'] ; % 0 NDF (BW 30-4-.48)
DataBlock(15).DsPath = [RootAnalysisDir,'2015-10-30-0/data007/data007'] ; % 0 NDF (2 TP, 1 SP)

% CNG KO resecue (TX at 2 m) from Jeannie Chen's lab
% missing triggers data000, 0001, ...
DataBlock(16).SpontPath{1} = [RootAnalysisDir,'2015-11-02-0/data000/data000'] ; % in darkness
DataBlock(16).SpontPath{2} = [RootAnalysisDir,'2015-11-02-0/data009/data009'] ; % ndf 0 grey screen
DataBlock(16).DfPath = [RootAnalysisDir,'2015-11-02-0/data001/data001'] ; % 3 NDF
DataBlock(16).DfParams.NDF =   [5,5,5,4,4,4,3,3,3,2,2,2,1,1,1] ; % on filter turret 
DataBlock(16).DfParams.Ftime = [2,4,8,2,4,8,2,4,8,2,4,8,2,4,8] ; % ms
DataBlock(16).DfParams.interFlashInt = [3] ; % sec
DataBlock(16).BwPath{1} = [RootAnalysisDir,'2015-11-02-0/data006/data006'] ; % 0 NDF (BW 10-1-.48)
DataBlock(16).BwPath{2} = [RootAnalysisDir,'2015-11-02-0/data007/data007'] ; % 0 NDF (BW 30-4-.48)
DataBlock(16).BwPath{3} = [RootAnalysisDir,'2015-11-02-0/data003/data003'] ; % 2 NDF (BW 10-2-.48 3600s)
DataBlock(16).BwPath{4} = [RootAnalysisDir,'2015-11-02-0/data004/data004'] ; % 2 NDF (BW 30-4-.48 3600s)
DataBlock(16).BwMoviePath{1} = [RootMovieDir,] ;
DataBlock(16).FfStepPath{1} = [RootAnalysisDir,'2015-11-02-0/data002/data002'] ; % 3 NDF gwgb
DataBlock(16).FfStepPath{2} = [RootAnalysisDir,'2015-11-02-0/data005/data005'] ; % 0 NDF gwgb
DataBlock(16).DsPath = [RootAnalysisDir,'2015-11-02-0/data008/data008'] ; % 0 NDF (2 TP, 2 SP, 12 directions)
DataBlock(16).mapEiFlag = 'true' ;

% CNG het chimera (~ 1 month old from Jeannie Chen's lab
% 80/20 NDF 3 but than switched to 50/50 NDF 0 
DataBlock(17).SpontPath{1} = [RootAnalysisDir,'2016-04-28-0/data000/data000'] ; % in darkness
DataBlock(17).DfPath = [RootAnalysisDir,'2016-04-28-0/data001-2/data001-2'] ; % 3 NDF
DataBlock(17).DfParams.NDF =   [5,5,5,5,4,4,4,4,4,4,5,3,3,3,3,3,3,3,2,2,2,1,1,1] ; % on filter turret
DataBlock(17).DfParams.Ftime = [2,4,8,4,2,4,8,2,4,8,8,2,4,8,2,4,8,4,2,4,8,2,4,8] ; % ms
DataBlock(17).DfParams.interFlashInt = 3 ; % sec
DataBlock(17).FfStepPath{1} = [RootAnalysisDir,'2016-04-28-0/data003/data003'] ; % 3 NDF gwgb
DataBlock(17).FfStepPath{2} = [RootAnalysisDir,'2016-04-28-0/data006/data006'] ; % 0 NDF gwgb
DataBlock(17).BwPath{1} = [RootAnalysisDir,'2016-04-28-0/data004/data004'] ; % 2 NDF (BW 10-2-.48)
DataBlock(17).BwPath{2} = [RootAnalysisDir,'2016-04-28-0/data005/data005'] ; % 2 NDF (BW 30-4-.48)
DataBlock(17).BwPath{3} = [RootAnalysisDir,'2016-04-28-0/data007/data007'] ; % 0 NDF (BW 10-1-.48)
DataBlock(17).BwPath{4} = [RootAnalysisDir,'2016-04-28-0/data008/data008'] ; % 0 NDF (BW 30-4-.48)
DataBlock(17).DsPath = [RootAnalysisDir,'2016-04-28-0/data009/data009'] ; % 0 NDF (2 TP, 1 SP)

% data with CNG experiment stimuli (dim flash, etc) on mice that were like
% wrong genotype (2016-04-27, 2016-04-28, 2016-05-13)

% C57/Cre-ERT2 control mouse for CNG-CreER
% 2016-05-10, 2016-05-18

% 2016-05-10 has issues - see notes

% Cre-ERT2/c57 Fed TM 
% naked retina
DataBlock(18).DfPath = [RootAnalysisDir,'2016-05-18-0/data000-1/data000-1'] ; % 3 NDF
DataBlock(18).DfParams.NDF =   [5,5,5,5,5,5,4,4,4,4,4,4,3,3,3,2,2,2,1,1,1] ; % on filter turret
DataBlock(18).DfParams.Ftime = [2,4,8,4,2,8,2,4,8,2,4,8,2,4,8,2,4,8,2,4,8] ; % ms
DataBlock(18).DfParams.interFlashInt = 3 ; % sec
DataBlock(18).FfStepPath{1} = [RootAnalysisDir,'2016-05-18-0/data002/data002'] ; % 3 NDF gwgb
DataBlock(18).FfStepPath{2} = [RootAnalysisDir,'2016-05-18-0/data003/data003'] ; % 2 NDF gwgb
DataBlock(18).FfStepPath{2} = [RootAnalysisDir,'2016-05-18-0/data006/data006'] ; % 2 NDF gwgb
DataBlock(18).BwPath{1} = [RootAnalysisDir,'2016-05-18-0/data004/data004'] ; % 2 NDF (BW 10-2-.48)
DataBlock(18).BwPath{2} = [RootAnalysisDir,'2016-05-18-0/data005/data005'] ; % 2 NDF (BW 30-4-.48)
DataBlock(18).BwPath{3} = [RootAnalysisDir,'2016-05-18-0/data007/data008'] ; % 0 NDF (BW 10-1-.48)
DataBlock(18).DsPath = [RootAnalysisDir,'2016-05-18-0/data007/data007'] ; % 0 NDF (2 TP, 1 SP)
DataBlock(18).DfBwBlock = 3 ; % BW block to use for classificiation 

% CNGB KO recorded ~3 months
% good recording
DataBlock(19).SpontPath{1} = [RootAnalysisDir,'2017-01-04-0/data000/data000'] ; % in darkness
DataBlock(19).DfPath = [RootAnalysisDir,'2017-01-04-0/data001_tw/data001_tw'] ; % 3 NDF
DataBlock(19).DfParams.NDF =   [5,5,5,4,4,4,4,3,3,3,3,3,3,2,2,2,2,2,2,1,1,1,1,1,1] ; % on filter turret
DataBlock(19).DfParams.Ftime = [8,4,8,8,4,2,8,8,4,2,4,2,8,2,4,8,2,4,8,8,2,4,2,4,8] ; % ms
DataBlock(19).DfParams.interFlashInt = 3 ; % sec
DataBlock(19).FfStepPath{1} = [RootAnalysisDir,'2017-01-04-0/data002/data002'] ; % 3 NDF gwgb
DataBlock(19).FfStepPath{2} = [RootAnalysisDir,'2017-01-04-0/data007/data007'] ; % 0 NDF gwgb
DataBlock(19).BwPath{1} = [RootAnalysisDir,'2017-01-04-0/data003/data003'] ; % 2 NDF (BW 10-2-.48) stopped short
DataBlock(19).BwPath{2} = [RootAnalysisDir,'2017-01-04-0/data004/data004'] ; % 2 NDF (BW 30-4-.48)
DataBlock(19).BwPath{3} = [RootAnalysisDir,'2017-01-04-0/data005_offset/data005_offset'] ; % 0 NDF (BW 15-2-.48)
DataBlock(19).DsPath = [RootAnalysisDir,'2017-01-04-0/data006/data006'] ; % 0 NDF (2 TP, 1 SP)
DataBlock(19).DfBwBlock = 3 ; % BW block to use for classificiation 
DataBlock(19).ExampleCelli = 15 ; % 258 (for nonspike sorted set) 

% CNGB KO recorded ~3 months
% good recording - but noise prevented a complete recording 
% good mosaics in BwPath{1}
DataBlock(20).SpontPath{1} = [RootAnalysisDir,'2017-01-05-0/data000/data000'] ; % in darkness
DataBlock(20).DfPath = [RootAnalysisDir,'2017-01-05-0/data001_tw/data001_tw'] ; % 3 NDF
DataBlock(20).DfParams.NDF =   [5,5,5,5,4,4,4,4,3,3,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1] ; % on filter turret
DataBlock(20).DfParams.Ftime = [8,2,4,8,8,2,4,8,8,2,4,8,4,2,2,4,2,4,8,2,4,8,2,4,8] ; % ms
DataBlock(20).DfParams.interFlashInt = 3 ; % sec
DataBlock(20).FfStepPath{1} = [RootAnalysisDir,'2017-01-05-0/data002/data002'] ; % 3 NDF gwgb
DataBlock(20).BwPath{1} = [RootAnalysisDir,'2017-01-05-0/data003/data003'] ; % 0 NDF (BW 10-1-.48) has significant noise contamination
DataBlock(20).DfBwBlock = 1 ; % BW block to use for classificiation 
DataBlock(20).ExampleCelli = 102 ;% 260 (non-spike sorted set)
DataBlock(20).ExampleFlashi = 10 ;

% CNGB KO recorded ~3 months
% ok recording but significant noise problems
% this retina had very weak light responses
DataBlock(21).SpontPath{1} = [RootAnalysisDir,'2017-01-06-0/data000/data000'] ; % in darkness
DataBlock(21).DfPath = [RootAnalysisDir,'2017-01-06-0/data001-2_tw/data001-2_tw'] ; % 3 NDF
DataBlock(21).DfParams.NDF =   [5,5,5,5,4,4,4,4,3,3,3,3,3,2,2,2,2,2,2,1,1,1,1,1,1] ; % on filter turret
DataBlock(21).DfParams.Ftime = [8,2,4,8,8,4,8,8,8,4,2,4,8,2,4,8,2,4,8,8,2,4,8,4,2] ; % ms
DataBlock(21).DfParams.interFlashInt = 3 ; % sec
DataBlock(21).FfStepPath{1} = [RootAnalysisDir,'2017-01-06-0/data003/data003'] ; % 3 NDF gwgb
DataBlock(21).FfStepPath{2} = [RootAnalysisDir,'2017-01-06-0/data006/data006'] ; % 0 NDF gwgb
DataBlock(21).BwPath{1} = [RootAnalysisDir,'2017-01-06-0/data004/data004'] ; % 2 NDF (BW 30-4-.48)
DataBlock(21).BwPath{2} = [RootAnalysisDir,'2017-01-06-0/data005/data005'] ; % 0 NDF (BW 30-4-.48)
DataBlock(21).DsPath = [RootAnalysisDir,'2017-01-06-0/data007/data007'] ; % 0 NDF (2 TP, 1 SP)
DataBlock(21).DfBwBlock = 2 ; % BW block to use for classificiation 
DataBlock(21).ExampleCelli = 1 ; % 45 (non-spike sorted set)

% CNGB -Cre negative fed TM @1m recorded ~3.5 months
% electrical noise issues - but data may still be useable
DataBlock(22).DfPath = [RootAnalysisDir,'2017-01-12-0/data000/data000'] ; % 3 NDF
DataBlock(22).DfParams.NDF =   [5,5,5,5,4,4,4,4,4,4,3,3,3,3,3,2,2,2,2,2,2,1,1,1,1,1] ; % on filter turret
DataBlock(22).DfParams.Ftime = [8,2,4,8,8,4,2,4,2,8,8,4,2,4,8,8,4,2,4,2,8,2,4,2,4,8] ; % ms
DataBlock(22).DfParams.interFlashInt = 3 ; % sec
DataBlock(22).FfStepPath{1} = [RootAnalysisDir,'2017-01-12-0/data001/data001'] ; % 3 NDF gwgb
DataBlock(22).BwPath{1} = [RootAnalysisDir,'2017-01-12-0/data002/data002'] ; % 0 NDF (BW 15-2-.48)
DataBlock(22).DsPath = [RootAnalysisDir,'2017-01-12-0/data003/data003'] ; % 0 NDF (2 TP, 1 SP)
DataBlock(22).DfBwBlock = 1 ; % BW block to use for classificiation 

% CNGB +Cre rescue fed TM @ 1m recorded ~3.5 months
DataBlock(23).SpontPath{1} = [RootAnalysisDir,'2017-01-19-0/data000/data000'] ; % in darkness
DataBlock(23).DfPath = [RootAnalysisDir,'2017-01-19-0/data001_JC/data001_JC'] ; % 3 NDF
DataBlock(23).DfParams.NDF =   [5,5,5,5,4,4,4,4,4,3,3,3,3,3,3,2,2,2,2,2,2,1,1,1] ; % on filter turret
DataBlock(23).DfParams.Ftime = [8,2,4,8,8,4,2,4,8,2,4,8,4,2,8,8,2,4,2,4,8,4,2,8] ; % ms
DataBlock(23).DfParams.interFlashInt = 3 ; % sec
DataBlock(23).FfStepPath{1} = [RootAnalysisDir,'2017-01-19-0/data002/data002'] ; % 3 NDF gwgb
DataBlock(23).BwPath{1} = [RootAnalysisDir,'2017-01-19-0/data003/data003'] ; % 0 NDF (BW 15-2-.48)
DataBlock(23).DfBwBlock = 1 ; % BW block to use for classificiation
DataBlock(23).ExampleCelli = 23 ; % 44 (non-spike sorted set)

% CNGB +Cre rescue fed TM @ 1m recorded ~3.5 months
DataBlock(24).SpontPath{1} = [RootAnalysisDir,'2017-01-20-0/data000/data000'] ; % in darkness (steadly better recording during this datablock)
DataBlock(24).DfPath = [RootAnalysisDir,'2017-01-20-0/data001_JC/data001_JC'] ; % 3 NDF
DataBlock(24).DfParams.NDF =   [5,5,5,5,5,4,4,4,4,4,4,3,3,3,3,3,3,2,2,2,2,2,2,1,1,1,1,1] ; % on filter turret
DataBlock(24).DfParams.Ftime = [2,4,8,4,8,8,2,4,8,4,2,2,4,2,8,4,8,8,2,4,8,2,4,4,2,2,8,4] ; % ms
DataBlock(24).DfParams.interFlashInt = 3 ; % sec
DataBlock(24).FfStepPath{1} = [RootAnalysisDir,'2017-01-20-0/data003/data003'] ; % 3 NDF gwgb
DataBlock(24).FfStepPath{2} = [RootAnalysisDir,'2017-01-20-0/data006/data006'] ; % 0 NDF gwgb
DataBlock(24).BwPath{1} = [RootAnalysisDir,'2017-01-20-0/data004/data004'] ; % 2 NDF (BW 30-4-.48)
DataBlock(24).BwPath{2} = [RootAnalysisDir,'2017-01-20-0/data005/data005'] ; % 0 NDF (BW 30-4-.48)
DataBlock(24).DsPath = [RootAnalysisDir,'2017-01-20-0/data007/data007'] ; % 0 NDF (2 TP, 1 SP)
DataBlock(24).DfBwBlock = 2 ; % BW block to use for classificiation 
DataBlock(24).ExampleCelli = 81 ;% 72 (non-spike sorted set)
DataBlock(24).ExampleFlashi = 10 ;

% CNGBdeltaCaM bred at Duke
% data003 has a sta delay
DataBlock(25).DfPath = [RootAnalysisDir,'2017-02-10-0/data000_JC/data000_JC'] ; % 3 NDF
DataBlock(25).DfParams.NDF =   [5,5,5,5,5,5,4,4,4,4,4,4,3,3,3,3,3,3,2,2,2,2,2,2,1,1,1,1,1] ; % on filter turret
DataBlock(25).DfParams.Ftime = [8,2,4,2,4,8,8,2,4,8,4,2,2,2,4,8,4,8,2,4,8,2,4,8,2,4,2,4,8] ; % ms
DataBlock(25).DfParams.interFlashInt = 3 ; % sec
DataBlock(25).FfStepPath{1} = [RootAnalysisDir,'2017-02-10-0/data001/data001'] ; % 3 NDF gwgb
DataBlock(25).BwPath{1} = [RootAnalysisDir,'2017-02-10-0/data002/data002'] ; % 2 NDF (BW 10-2-.48)
DataBlock(25).BwPath{2} = [RootAnalysisDir,'2017-02-10-0/data003/data003'] ; % 0 NDF (BW 10-1-.48)
DataBlock(25).DfBwBlock = 1 ; % BW block to use for classificiation 
DataBlock(25).ExampleCelli = 28 ; %11 (non-spike sorted set)
DataBlock(25).ExampleFalshi = 10 ;

% CNGB 1m KO bred at Duke
% No ON cells in STA at rod light levels, no STA structure at cone light level
DataBlock(26).DfPath = [RootAnalysisDir,'2017-02-13-0/data000_JC/data000_JC'] ; % 3 NDF
DataBlock(26).DfParams.NDF =   [5,5,5,5,4,4,4,4,4,4,3,3,3,3,3,2,2,2,2,2,1,1,1,1,1] ; % on filter turret
DataBlock(26).DfParams.Ftime = [8,2,4,8,8,2,4,8,2,4,4,2,8,4,8,2,4,8,4,2,2,4,2,4,8] ; % ms
DataBlock(26).DfParams.interFlashInt = 3 ; % sec
DataBlock(26).FfStepPath{1} = [RootAnalysisDir,'2017-02-13-0/data001/data001'] ; % 3 NDF gwgb
DataBlock(26).BwPath{1} = [RootAnalysisDir,'2017-02-13-0/data002/data002'] ; % 2 NDF (BW 10-2-.48)
DataBlock(26).BwPath{2} = [RootAnalysisDir,'2017-02-13-0/data003/data003'] ; % 0 NDF (BW 10-1-.48)
DataBlock(26).DsPath = [RootAnalysisDir,'2017-02-13-0/data005/data005'] ; % 0 NDF (2 TP, 1 SP)
DataBlock(26).DfBwBlock = 1 ; % BW block to use for classificiation 
DataBlock(26).ExampleCelli = 1 ; % random 

% CNGBdeltaCaM bred at Duke
DataBlock(27).DfPath = [RootAnalysisDir,'2017-02-14-0/data000-1_JC/data000-1_JC'] ; % 3 NDF
DataBlock(27).DfParams.NDF =   [5,5,5,5,5,5,4,4,4,4,4,4,4,3,3,3,3,3,3,3,2,2,2,2,2,2,1,1,1,1,1] ; % on filter turret
DataBlock(27).DfParams.Ftime = [8,4,2,4,8,2,2,4,8,4,2,8,2,2,4,8,2,2,4,8,8,2,4,2,8,8,2,4,2,4,8] ; % ms
DataBlock(27).DfParams.interFlashInt = 3 ; % sec
DataBlock(27).BwPath{1} = [RootAnalysisDir,'2017-02-14-0/data002/data002'] ; % 2 NDF (BW 10-1-.48)
DataBlock(27).DfBwBlock = 1 ; % BW block to use for classificiation 
DataBlock(27).ExampleCelli = 77 ; % random

% CNGB 1m KO bred at Duke
DataBlock(28).DfPath = [RootAnalysisDir,'2017-02-16-0/data000_JC/data000_JC'] ; % 3 NDF
DataBlock(28).DfParams.NDF =   [5,5,5,5,5,4,4,4,4,4,4,3,3,3,3,3,3,2,2,2,2,2,2,1,1,1,1,1,1] ; % on filter turret
DataBlock(28).DfParams.Ftime = [8,4,2,4,8,8,2,4,8,2,4,4,8,2,4,8,2,2,4,8,2,4,8,2,4,2,8,4,8] ; % ms
DataBlock(28).DfParams.interFlashInt = 3 ; % sec
DataBlock(28).FfStepPath{1} = [RootAnalysisDir,'2017-02-16-0/data001/data001'] ; % 3 NDF gwgb
DataBlock(28).BwPath{1} = [RootAnalysisDir,'2017-02-16-0/data002/data002'] ; % 0 NDF (BW 10-1-.48)
DataBlock(28).DsPath = [RootAnalysisDir,'2017-02-16-0/data003/data003'] ; % 0 NDF (2 TP, 1 SP)
DataBlock(28).DfBwBlock = 1 ; % BW block to use for classificiation 
DataBlock(28).ExampleCelli = 94 ; % 84 (non-spike sorted set)
DataBlock(28).ExampleFlashi = 10 ;

% CNGB 1m KO bred at Duke
DataBlock(29).DfPath = [RootAnalysisDir,'2017-02-17-0/data000_JC/data000_JC'] ; % 3 NDF
DataBlock(29).DfParams.NDF =   [5,5,5,5,5,5,4,4,4,4,4,4,3,3,3,3,3,3,2,2,2,2,2,2,1,1,1,1] ; % on filter turret
DataBlock(29).DfParams.Ftime = [8,2,4,8,4,2,2,4,8,2,4,8,8,2,4,2,4,8,2,4,8,2,8,4,4,2,4,8] ; % ms
DataBlock(29).DfParams.interFlashInt = 3 ; % sec
DataBlock(29).BwPath{1} = [RootAnalysisDir,'2017-02-17-0/data001/data001'] ; % 2 NDF (BW 10-1-.48)
DataBlock(29).DfBwBlock = 1 ; % BW block to use for classificiation 
DataBlock(29).ExampleCelli = 1 ; % random

% CNGB 3m KO bred at Duke
DataBlock(30).SpontPath{1} = [RootAnalysisDir,'2017-04-13-0/data000/data000'] ; % in darkness
DataBlock(30).SpontPath{2} = [RootAnalysisDir,'2017-04-13-0/data004/data004'] ; % NDF bg=.5 (~5000Rh*/sec)
DataBlock(30).DfPath = [RootAnalysisDir,'2017-04-13-0/data001_JC/data001_JC'] ; % 3 NDF
DataBlock(30).DfParams.NDF =   [5,5,5,5,5,4,4,4,4,4,4,3,3,3,3,3,3,2,2,2,2,2,2,1,1,1,1,1,1] ; % on filter turret
DataBlock(30).DfParams.Ftime = [8,4,2,4,8,8,4,2,4,8,2,2,4,8,4,8,2,2,4,8,2,4,8,2,4,8,4,8,2] ; % ms
DataBlock(30).DfParams.interFlashInt = 3 ; % sec
DataBlock(30).BwPath{1} = [RootAnalysisDir,'2017-04-13-0/data002/data002'] ; % 0 NDF (BW 10-1-.48)
DataBlock(30).DfBwBlock = 1 ; % BW block to use for classificiation 
DataBlock(30).ExampleCelli = 1 ; % random
DataBlock(30).DsPath = [RootAnalysisDir,'2017-04-13-0/data003/data003'] ; % 0 NDF (2 TP, 1 SP)

% CNGB 3m KO bred at Duke
DataBlock(31).SpontPath{1} = [RootAnalysisDir,'2017-04-18-0/data000/data000'] ; % in darkness
DataBlock(31).SpontPath{2} = [RootAnalysisDir,'2017-04-18-0/data004/data004'] ; % NDF bg=.5 (~5000Rh*/sec)
DataBlock(31).DfPath = [RootAnalysisDir,'2017-04-18-0/data001_JC/data001_JC'] ; % 3 NDF
DataBlock(31).DfParams.NDF =   [5,5,5,5,5,5,4,4,4,4,4,4,3,3,3,3,3,3,2,2,2,2,2,2,1,1,1,1,1] ; % on filter turret
DataBlock(31).DfParams.Ftime = [2,4,8,4,2,8,8,2,4,8,4,2,2,4,8,2,4,8,8,4,2,8,4,2,2,4,8,2,4] ; % ms
DataBlock(31).DfParams.interFlashInt = 3 ; % sec
DataBlock(31).BwPath{1} = [RootAnalysisDir,'2017-04-18-0/data002/data002'] ; % 0 NDF (BW 10-1-.48)
DataBlock(31).DfBwBlock = 1 ; % BW block to use for classificiation 
DataBlock(31).ExampleCelli = 66 ; % 
DataBlock(31).DsPath = [RootAnalysisDir,'2017-04-18-0/data003/data003'] ; % 0 NDF (2 TP, 1 SP)

% ELFN1 from Kirill M Lab (RRL mark)
% sevaral noise events in the dim flash and spont data
% generally good recording
DataBlock(32).SpontPath{1} = [RootAnalysisDir,'2017-05-30-0/data000/data000'] ; % in darkness
DataBlock(32).DfPath = [RootAnalysisDir,'2017-05-30-0/data001/data001'] ; % 3 NDF
DataBlock(32).DfParams.NDF =   [5,5,5,5,5,4,4,4,4,4,4,3,3,3,3,3,3,2,2,2,2,2,2,1,1,1,1,1,1] ; % on filter turret
DataBlock(32).DfParams.Ftime = [8,4,2,8,4,4,2,8,4,2,8,8,4,2,4,8,2,2,4,8,2,4,8,8,2,4,2,8,4] ; % ms
DataBlock(32).DfParams.interFlashInt = 3 ; % sec
DataBlock(32).BwPath{1} = [RootAnalysisDir,'2017-05-30-0/data005/data005'] ; % 2 NDF (BW 30-4-.48)
DataBlock(32).BwPath{2} = [RootAnalysisDir,'2017-05-30-0/data006/data006'] ; % 0 NDF (BW 10-1-.48)
DataBlock(32).DfBwBlock = 2 ; % BW block to use for classificiation 
DataBlock(32).ExampleCelli = 1 ; % 
DataBlock(32).DsPath = [RootAnalysisDir,'2017-05-30-0/data007/data007'] ; % 0 NDF (2 TP, 1 SP)
DataBlock(32).FfStepPath{1} = [RootAnalysisDir,'2017-05-30-0/data003/data003'] ; % 3 NDF gwgb

% ELFN1 from Kirill M Lab (cage KM 647 L ear mark)
% not a great recording
DataBlock(33).DfPath = [RootAnalysisDir,'2017-06-28-0/data000/data000'] ; % 3 NDF
DataBlock(33).DfParams.NDF =   [5,4,3,2,2,2,2,2,2,1,1,1,1,1,1,0,0,0] ; % on filter turret
DataBlock(33).DfParams.Ftime = [8,8,8,8,2,4,2,8,4,2,4,2,8,4,8,2,4,2] ; % ms
DataBlock(33).DfParams.interFlashInt = 3 ; % sec
DataBlock(33).BwPath{1} = [RootAnalysisDir,'2017-06-28-0/data001/data001'] ; % 0 NDF (BW 15-2)
DataBlock(33).DfBwBlock = 1 ; % BW block to use for classificiation 
DataBlock(33).ExampleCelli = 1 ; % 
DataBlock(33).DsPath = [RootAnalysisDir,'2017-06-28-0/data002/data002'] ; % 0 NDF (2 TP, 1 SP)

% ELFN1 from Kirill M Lab (cage 647 RR mark)
% generally good recording
DataBlock(34).SpontPath{1} = [RootAnalysisDir,'2017-07-12-0/data000-1/data000-1'] ; % in darkness
DataBlock(34).DfPath = [RootAnalysisDir,'2017-07-12-0/data002/data002'] ; % 3 NDF
DataBlock(34).DfParams.NDF =   [5,4,3,3,3,3,2,2,2,2,2,2,1,1,1,1,1,1,0] ; % on filter turret
DataBlock(34).DfParams.Ftime = [8,8,8,4,8,4,2,4,2,8,4,8,2,4,8,2,4,8,2] ; % ms
DataBlock(34).DfParams.interFlashInt = 3 ; % sec
DataBlock(34).BwPath{1} = [RootAnalysisDir,'2017-07-12-0/data003/data003'] ; % 2 NDF (BW 15-2-.48)
DataBlock(34).BwPath{2} = [RootAnalysisDir,'2017-07-12-0/data004/data004'] ; % 2 NDF (BW 30-4-.48)
DataBlock(34).BwPath{3} = [RootAnalysisDir,'2017-07-12-0/data005/data005'] ; % 0 NDF (BW 10-1-.48) - no ON cells
DataBlock(34).DfBwBlock = 2 ; % BW block to use for classificiation 
DataBlock(34).ExampleCelli = 1 ; % 
DataBlock(34).DsPath = [RootAnalysisDir,'2017-07-12-0/data006/data006'] ; % 0 NDF (2 TP, 1 SP)
DataBlock(34).RgPath = [RootAnalysisDir,'2017-07-12-0/data007/data007'] ; % 0 NDF reverse grating 50% contrast 3 phases, 1 TP, 6 SP
DataBlock(34).RgPath = [RootAnalysisDir,'2017-07-12-0/data008/data008'] ; % 0 NDF reverse grating 100% contrast 3 phases, 1 TP, 6 SP
DataBlock(34).RepStim{1} = [RootAnalysisDir,'2017-07-12-0/data009/data009'] ; % repeated cat movie


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


function [DataBlock,Params] = DataBlocks_NaturalStim 

% DataBlocks for natural stimulation movies (run as script or load as function)

%% functions to perform
perform.MovieRasterPlots = false ; 

perform.MovieSpikeParamFinder = false;

perform.DsCellFinder = false ;

perform.DsAnalysisNatStim = false ;
perform.DsAnalysisNatStimV2 = false ;
perform.DsAnalysisNatStimV3 = false ;

perform.TformCalculator = false ; % calculate a Tform for EI-->monitor coordinates (not an inpendent function)

perform.SquareAnalysisV2 = false ;

perform.NatStimMappedBw = false ;

perform.ImJitterAnalysis = false ;
perform.ImJitterAnalysisV2 = false ;

perform.RepStimNecSufAnalysis = false ;
perform.RepStimNecSufAnalysisV2 = false ;

perform.ImSliderAnalysis = false ;

perform.ImSliderConcatAnalysis = false ;
perform.ImSliderConcatAnalysisV2 = false ;
perform.ImSliderConcatAnalysisV3 = false ;
perform.ImSliderConcatAnalysisV4 = false ;
perform.ImSliderConcatAnalysisV5 = false ;
perform.ImSliderConcatAnalysisV6 = false ;

perform.RepStimClusterAnalysis = false ;

perform.RevMovConcatAnalyzer = false ;

%% perform functions on this data block set
DBset_id = 1 ;

% data block sets
DBset{1} = [19] ; % temp analysis
DBset{2} = [1:4] ; % Cells with...
Params.MovieNum = 1 ; % {index} of movie to analyze
Params.DsPathNum = 1 ;
ForIgor.exportFlag = 0 ; % can't export logicals

%% Data Blocks defined
% Root dircetory for analysis files
%RootAnalysisDir = '/Volumes/Berlioz/Analysis/' ; % Berlioz
RootAnalysisDir = '/Volumes/lab/Experiments/Array/Analysis/' ; % Brahms server
RootMovieDir = '/Volumes/lab/acquisition/movie-xml/' ; % Brahms server xml movie root
RootImageDir = '/Volumes/lab/Experiments/Array/Images/' ; 
RootSavePath = '/Users/jcafaro/Documents/AnalysisFigures/' ;
RootNatStimMoviesDir = '/Volumes/lab/Documents/Movies/' ; % natural stim moves

% c57, sparse array, not very sensitive mount
% cannot identify ds cells in dg
DataBlock(1).BwPath{1} = [RootAnalysisDir,'2016-10-14-0/data000/data000'] ; % 
DataBlock(1).MovieRepPath{1} = [RootAnalysisDir,'2016-10-14-0/data002/data002'] ; % mouse movie
DataBlock(1).MovieRepPath{2} = [RootAnalysisDir,'2016-10-14-0/data003/data003'] ; % squirrel movie
DataBlock(1).MovieRepPath{3} = [RootAnalysisDir,'2016-10-14-0/data004/data004'] ; % cat movie
DataBlock(1).MovieRepPath{4} = [RootAnalysisDir,'2016-10-14-0/data005/data005'] ; % cat movie slow
DataBlock(1).MovieRepPath{5} = [RootAnalysisDir,'2016-10-14-0/data006/data006'] ; % cat movie zoom
DataBlock(1).MovieRepPath{6} = [RootAnalysisDir,'2016-10-14-0/data007/data007'] ; % cat movie phase scrambled
DataBlock(1).MovieRepPath{7} = [RootAnalysisDir,'2016-10-14-0/data008/data008'] ; % cat movie temp scrambled 
DataBlock(1).MovieRepPath{8} = [RootAnalysisDir,'2016-10-14-0/data009/data009'] ; % mouse movie again
DataBlock(1).MovieRepPath{9} = [RootAnalysisDir,'2016-10-14-0/data010/data010'] ; % cat movie spatial scramble
DataBlock(1).MovieRepPath{10} = [RootAnalysisDir,'2016-10-14-0/data011/data011'] ; % cat movie slow and zoom
DataBlock(1).MovieRepFrameNum = [1300*4, 2586*2, 4999, 1245*4, 4999, 4999, 4999, 1300*4, 4999, 4999]  ; 
DataBlock(1).FfPulse = [RootAnalysisDir,'2016-10-14-0/data002/data002'] ; %
DataBlock(1).DsPath{1} = [RootAnalysisDir,'2016-10-14-0/data001/data001'] ; % low mean (2 TP, 8 directions)
DataBlock(1).DsPath{2} = [RootAnalysisDir,'2016-10-14-0/data012/data012'] ; % higher mean (2 TP, 8 directions)

% c57, dense array, ok mount
DataBlock(2).BwPath{1} = [RootAnalysisDir,'2016-10-26-0/data000/data000'] ; % 15-2 53x40
DataBlock(2).MovingSquare = [RootAnalysisDir,'2016-10-26-0/data001/data001'] ; % 15 width, full shift
DataBlock(2).MovieDsConcat{1} = [RootAnalysisDir,'2016-10-26-0/data002-3/data002-3'] ;% squirrel
DataBlock(2).MovieDsConcat{2} = [RootAnalysisDir,'2016-10-26-0/data002-6/data002-6'] ; % cat
DataBlock(2).MovieDsConcat{3} = [RootAnalysisDir,'2016-10-26-0/data002-7/data002-7'] ; % mouse
DataBlock(2).MovieRepPath{1} = [RootAnalysisDir,'2016-10-26-0/data003/data003'] ; % full squirrel movie 15x
DataBlock(2).MovieRepPath{2} = [RootAnalysisDir,'2016-10-26-0/data006/data006'] ; % full cat movie 15x
DataBlock(2).MovieRepPath{3} = [RootAnalysisDir,'2016-10-26-0/data007/data007'] ; % full mouse movie 15x
DataBlock(2).MovieRepPath{4} = [RootAnalysisDir,'2016-10-26-0/data005/data005'] ; % squirrel movie (10s 100x)
DataBlock(2).MovieRepPath{5} = [RootAnalysisDir,'2016-10-26-0/data009/data009'] ; % cat movie (10s 100x)
DataBlock(2).MovieRepFrameNum = [2586*2, 4999, 1300*4] ; % number of frames for each movie
DataBlock(2).MovieStixWidth = [1,2,1] ; % stixel width
DataBlock(2).MovieFrameInterval = [2,1,4] ; % number of frames for each movie
DataBlock(2).xstart = [0,80,80] ;
DataBlock(2).ystart = [0,60,60] ;
DataBlock(2).BwRepeat{1} = [RootAnalysisDir,'2016-10-26-0/data008/data008'] ; % 15-2 53x40
DataBlock(2).GwnFf = [RootAnalysisDir,'2016-10-26-0/data004/data004'] ; % large stixel gaussian white noise (effectively full field)
DataBlock(2).DsPath{1} = [RootAnalysisDir,'2016-10-26-0/data002/data002'] ; % (2 TP, 8 directions)
DataBlock(2).BwImagePath = [RootImageDir,'2016-10-26-0/4x OLED mapping 1.jpg'] ; % BW 15-2 53x40
DataBlock(2).ArrayImagePath = [RootImageDir,'2016-10-26-0/4x piece IR.jpg'] ; % 
DataBlock(2).BwMoviePath = [RootMovieDir,'BW-15-2-0.48-11111-53x40-60.35.xml'] ; %
DataBlock(2).TformEiPath = [RootSavePath,'Db2TformEi'] ; % path of saved TformEi
DataBlock(2).MoviePath{1} = [RootNatStimMoviesDir,'squirrel_video/squirrel_mean117_sd62_0to255.mat'] ; % squirrel movie
DataBlock(2).MoviePath{2} = [RootNatStimMoviesDir,'CatCam/cat_mean117_sd62_0to255.mat'] ; % cat movie
DataBlock(2).MoviePath{3} = [RootNatStimMoviesDir,'MrscFlogelMouseCam/mouse_mean117_sd62_0to255.mat'] ; % cat movie

% c57, dense array, ok mount
DataBlock(3).BwPath{1} = [RootAnalysisDir,'2016-11-18-0/data000/data000'] ; % 15-2 53x40
DataBlock(3).MovingSquare = [RootAnalysisDir,'2016-11-18-0/data001/data001'] ; % 15 width, full shift
DataBlock(3).MovieDsConcat{1} = [RootAnalysisDir,'2016-11-18-0/data002-3/data002-3'] ;% moving bar movie fast
DataBlock(3).MovieDsConcat{2} = [RootAnalysisDir,'2016-11-18-0/data002-4/data002-4'] ; % moving bar movie slow
DataBlock(3).MovieDsConcat{3} = [RootAnalysisDir,'2016-11-18-0/data002-5/data002-5'] ; % squirrel movie
DataBlock(3).MovieDsConcat{4} = [RootAnalysisDir,'2016-11-18-0/data006-8/data006-8'] ; % cat movie
DataBlock(3).MovieDsConcat{5} = [RootAnalysisDir,'2016-11-18-0/data006-9/data006-9'] ; % moving bar movie high contrast
DataBlock(3).MovieRepFrameNum = [1959,1959*2,2586*2,4999,1959*2] ; % number of frames for each movie
DataBlock(3).MovieStixWidth = [2,2,1,2,2] ; % stixel width
DataBlock(3).MovieFrameInterval = [1,2,2,1,2] ; % number of frames for each movie
DataBlock(3).xstart = [0,0,0,80,0] ;
DataBlock(3).ystart = [0,0,0,60,0] ;
DataBlock(3).DsPath{1} = [RootAnalysisDir,'2016-11-18-0/data002/data002'] ; % (2 TP, 8 directions) - use for Movies 1-3
DataBlock(3).DsPath{2} = [RootAnalysisDir,'2016-11-18-0/data006/data006'] ; %(2 TP, 8 directions)- use for Movies 4-5    
DataBlock(3).BwImagePath = [RootImageDir,'2016-11-18-0/4x OLED 1.jpg'] ; % BW 15-2 53x40
DataBlock(3).ArrayImagePath = [RootImageDir,'2016-11-18-0/4x IR.jpg'] ; % 
DataBlock(3).BwMoviePath = [RootMovieDir,'BW-50-2-0.48-11111-16x12.xml'] ; %
DataBlock(3).TformEiPath = [RootSavePath,'Db3TformEi'] ; % path of saved TformEi
DataBlock(3).MoviePath{1} = [RootNatStimMoviesDir,'moving_bars/MovingBarsMovSqueezed.mat'] ; % moving bar movie
DataBlock(3).MoviePath{2} = [RootNatStimMoviesDir,'moving_bars/MovingBarsMovSqueezed.mat'] ; % moving bar movie
DataBlock(3).MoviePath{3} = [RootNatStimMoviesDir,'squirrel_video/squirrel_mean117_sd62_0to255.mat'] ; % squirrel movie
DataBlock(3).MoviePath{4} = [RootNatStimMoviesDir,'CatCam/cat_mean117_sd62_0to255.mat'] ; % cat movie
DataBlock(3).MoviePath{5} = [RootNatStimMoviesDir,'moving_bars/MovingBarsMovHighContrast.mat'] ; % high contrast moving bar movie

% c57, dense array, 
% good mount early but lost areas some later 
DataBlock(4).BwPath{1} = [RootAnalysisDir,'2016-12-07-0/data000/data000'] ; % 20-1 40x30
DataBlock(4).MovieDsConcat{1} = [RootAnalysisDir,'2016-12-07-0/data001-2/data001-2'] ; % local motion movie :bobCat 
DataBlock(4).MovieDsConcat{2} = [RootAnalysisDir,'2016-12-07-0/data001-3/data001-3'] ; % local motion movie :Rats 
DataBlock(4).MovieDsConcat{3} = [RootAnalysisDir,'2016-12-07-0/data001-4/data001-4'] ; % local motion movie :CheetahApproach 
DataBlock(4).MovieDsConcat{4} = [RootAnalysisDir,'2016-12-07-0/data001-5/data001-5'] ; % local motion movie :Lions
DataBlock(4).MovieDsConcat{5} = [RootAnalysisDir,'2016-12-07-0/data001-6/data001-6'] ; % local motion movie :Bears
DataBlock(4).MovieDsConcat{6} = [RootAnalysisDir,'2016-12-07-0/data001-7/data001-7'] ; % local motion movie :CheetahFollow
DataBlock(4).MovieDsConcat{7} = [RootAnalysisDir,'2016-12-07-0/data001-8/data001-8'] ; % local motion movie :Mouse
DataBlock(4).MovieDsConcat{8} = [RootAnalysisDir,'2016-12-07-0/data001-9/data001-9'] ; % local motion movie :Rats smaller
DataBlock(4).MovieDsConcat{9} = [RootAnalysisDir,'2016-12-07-0/data010-2/data010-2'] ; % local motion movie :bobCat 
DataBlock(4).MovieDsConcat{10} = [RootAnalysisDir,'2016-12-07-0/data010-3/data010-3'] ; % local motion movie :Rats 
DataBlock(4).MovieDsConcat{11} = [RootAnalysisDir,'2016-12-07-0/data010-4/data010-4'] ; % local motion movie :CheetahApproach 
DataBlock(4).MovieDsConcat{12} = [RootAnalysisDir,'2016-12-07-0/data010-5/data010-5'] ; % local motion movie :Lions
DataBlock(4).MovieDsConcat{13} = [RootAnalysisDir,'2016-12-07-0/data010-6/data010-6'] ; % local motion movie :Bears
DataBlock(4).MovieDsConcat{14} = [RootAnalysisDir,'2016-12-07-0/data010-7/data010-7'] ; % local motion movie :CheetahFollow
DataBlock(4).MovieDsConcat{15} = [RootAnalysisDir,'2016-12-07-0/data010-8/data010-8'] ; % local motion movie :Mouse
DataBlock(4).MovieDsConcat{16} = [RootAnalysisDir,'2016-12-07-0/data010-9/data010-9'] ; % local motion movie :Rats smaller
DataBlock(4).MovieDsConcat{17} = [RootAnalysisDir,'2016-12-07-0/data010-11/data010-11'] ; % moving bar high contrast
DataBlock(4).MovieDsConcat{18} = [RootAnalysisDir,'2016-12-07-0/data010-12/data010-12'] ; % cat cam
DataBlock(4).MovieDsConcat{19} = [RootAnalysisDir,'2016-12-07-0/data010-13/data010-13'] ; % mouse cam
DataBlock(4).MovieDsConcat{20} = [RootAnalysisDir,'2016-12-07-0/data010-14/data010-14'] ; % squirrel cam
DataBlock(4).MovieRepFrameNum = [664*2,450*2,263,239*2,471*2,625*2,114*2,450*2,664*2,450*2,263,239*2,471*2,625*2,114*2,450*2,1959*2,4999,1300*4,2586*2] ; % number of frames for each movie
DataBlock(4).MovieStixWidth = [1,2,2,2,1,1,1,1,1,2,2,2,1,1,1,1,2,2,1,1] ; % stixel width
DataBlock(4).MovieFrameInterval = [2,2,1,2,2,2,2,2,2,2,1,2,2,2,2,2,2,1,4,2] ; % number of frames for each movie
DataBlock(4).xstart = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,80,0,0] ;
DataBlock(4).ystart = [100,0,0,0,100,100,100,100,100,0,0,0,100,100,100,100,0,60,0,0] ;
DataBlock(4).DsPath{1} = [RootAnalysisDir,'2016-12-07-0/data001/data001'] ; % (2 TP, 8 directions) - use for Movies 1-9
DataBlock(4).DsPath{2} = [RootAnalysisDir,'2016-12-07-0/data010/data010'] ; %(2 TP, 8 directions)- use for Movies 9-20    
DataBlock(4).BwImagePath = [RootImageDir,'2016-12-07-0/4x OLED image 1.jpg'] ; % BW 50
DataBlock(4).ArrayImagePath = [RootImageDir,'2016-12-07-0/4x IR image 1.jpg'] ; % 
DataBlock(4).BwMoviePath = [RootMovieDir,'BW-50-2-0.48-11111-16x12.xml'] ; %
DataBlock(4).TformEiPath = [RootSavePath,'Db4TformEi'] ; % path of saved TformEi
DataBlock(4).MoviePath{1} = [RootNatStimMoviesDir,'localMotionClips/BobCat.mat'] ; %
DataBlock(4).MoviePath{2} = [RootNatStimMoviesDir,'localMotionClips/Rats.mat'] ; %
DataBlock(4).MoviePath{3} = [RootNatStimMoviesDir,'localMotionClips/CheetahApproach.mat'] ; %
DataBlock(4).MoviePath{4} = [RootNatStimMoviesDir,'localMotionClips/Lions.mat'] ; %
DataBlock(4).MoviePath{5} = [RootNatStimMoviesDir,'localMotionClips/Bears.mat'] ; %
DataBlock(4).MoviePath{6} = [RootNatStimMoviesDir,'localMotionClips/CheetahFollow.mat'] ; %
DataBlock(4).MoviePath{7} = [RootNatStimMoviesDir,'localMotionClips/Mouse.mat'] ; %
DataBlock(4).MoviePath{8} = [RootNatStimMoviesDir,'localMotionClips/Rats.mat'] ; %
DataBlock(4).MoviePath{9} = [RootNatStimMoviesDir,'localMotionClips/BobCat.mat'] ; %
DataBlock(4).MoviePath{10} = [RootNatStimMoviesDir,'localMotionClips/Rats.mat'] ; %
DataBlock(4).MoviePath{11} = [RootNatStimMoviesDir,'localMotionClips/CheetahApproach.mat'] ; %
DataBlock(4).MoviePath{12} = [RootNatStimMoviesDir,'localMotionClips/Lions.mat'] ; %
DataBlock(4).MoviePath{13} = [RootNatStimMoviesDir,'localMotionClips/Bears.mat'] ; %
DataBlock(4).MoviePath{14} = [RootNatStimMoviesDir,'localMotionClips/CheetahFollow.mat'] ; %
DataBlock(4).MoviePath{15} = [RootNatStimMoviesDir,'localMotionClips/Mouse.mat'] ; %
DataBlock(4).MoviePath{16} = [RootNatStimMoviesDir,'localMotionClips/Rats.mat'] ; %
DataBlock(4).MoviePath{17} = [RootNatStimMoviesDir,'moving_bars/MovingBarsMovHighContrast.mat'] ; % high contrast moving bar movie
DataBlock(4).MoviePath{18} = [RootNatStimMoviesDir,'CatCam/cat_mean117_sd62_0to255.mat'] ; % cat cam
DataBlock(4).MoviePath{19} = [RootNatStimMoviesDir,'MrscFlogelMouseCam/mouse_mean117_sd62_0to255.mat'] ; % mouse cam
DataBlock(4).MoviePath{20} = [RootNatStimMoviesDir,'squirrel_video/squirrel_mean117_sd62_0to255.mat'] ; % squirrel cam

% cx36 KO
% XY experiment 
% NOT CORRECT TformEiPath
DataBlock(5).MovieDsConcat{1} = [RootAnalysisDir,'2016-12-12-0/data011-010/data011-010'] ; % local motion movie :bobCat 
DataBlock(5).MovieRepFrameNum = [2586*2] ; % number of frames for each movie
DataBlock(5).MovieStixWidth = [1] ; % stixel width
DataBlock(5).MovieFrameInterval = [2] ; % number of frames for each movie
DataBlock(5).DsPath{1} = [RootAnalysisDir,'2016-12-12-0/data011/data011'] ; %(2 TP, 8 directions)- use for Movies 9-20    
DataBlock(5).BwImagePath = [RootImageDir,'2016-12-07-0/4x OLED image 1.jpg'] ; % BW 50
DataBlock(5).ArrayImagePath = [RootImageDir,'2016-12-07-0/4x IR image 1.jpg'] ; % 
DataBlock(5).TformEiPath = [RootSavePath,'Db4TformEi'] ; % path of saved TformEi
DataBlock(5).MoviePath{1} = [RootNatStimMoviesDir,'squirrel_video/squirrel_mean117_sd62_0to255.mat'] ; % squirrel cam

% c57 image jitter
% all blocks are short (<10 minutes) because of memmorry issues
DataBlock(6).DsPath{1} = [RootAnalysisDir,'2016-12-20-0/data000/data000'] ;
DataBlock(6).ImJitterConcat{1} = [RootAnalysisDir,'2016-12-20-0/data000-2/data000-2'] ; %image (trees and building)
DataBlock(6).ImJitterConcat{2} = [RootAnalysisDir,'2016-12-20-0/data000-3/data000-3'] ; %image (trees and building)
DataBlock(6).ImJitterConcat{3} = [RootAnalysisDir,'2016-12-20-0/data000-4/data000-4'] ; %square (bg trees and building) 
DataBlock(6).ImJitterConcat{4} = [RootAnalysisDir,'2016-12-20-0/data006-5/data006-5'] ; %image (forest bright)
DataBlock(6).ImJitterConcat{5} = [RootAnalysisDir,'2016-12-20-0/data006-7/data006-7'] ; %square (bg forest bright)
DataBlock(6).ImJitterConcat{6} = [RootAnalysisDir,'2016-12-20-0/data006-8/data006-8'] ; %image (hole in sand)

% c57 image jitter
% very poor ds selection
% memorry issue is fixed but frame rate is NOT ~60 hz
DataBlock(7).DsPath{1} = [RootAnalysisDir,'2017-02-01-0/data000/data000'] ;
DataBlock(7).ImJitterConcat{1} = [RootAnalysisDir,'2017-02-01-0/data000-1/data000-1'] ; %image (trees and buildings)
DataBlock(7).ImJitterConcat{2} = [RootAnalysisDir,'2017-02-01-0/data000-2/data000-2'] ; %square

% c57 image jitter
% memorry issue is fixed but frame rate is NOT ~60 hz
% this response sesetivity is better in later trials
DataBlock(8).DsPath{2} = [RootAnalysisDir,'2017-03-08-0/data002/data002'] ;
DataBlock(8).DsPath{1} = [RootAnalysisDir,'2017-03-08-0/data005/data005'] ;
DataBlock(8).ImJitterConcat{3} = [RootAnalysisDir,'2017-03-08-0/data002-1/data002-1'] ; %image (trees and buildings)
DataBlock(8).ImJitterConcat{2} = [RootAnalysisDir,'2017-03-08-0/data002-3/data002-3'] ; %square (no bg - trees and builings)
DataBlock(8).ImJitterConcat{1} = [RootAnalysisDir,'2017-03-08-0/data005-4_TW/data005-4_TW'] ; %image (bright forest)
DataBlock(8).ImJitterConcat{4} = [RootAnalysisDir,'2017-03-08-0/data005-6/data005-6'] ; %square (no bg - bright forest)

% c57 image jitter
% noise issues
DataBlock(9).DsPath{1} = [RootAnalysisDir,'2017-03-17-0/data000/data000'] ;
DataBlock(9).DsPath{2} = [RootAnalysisDir,'2017-03-17-0/data003/data003'] ;
DataBlock(9).DsPath{3} = [RootAnalysisDir,'2017-03-17-0/data004/data004'] ;
DataBlock(9).ImJitterConcat{1} = [RootAnalysisDir,'2017-03-17-0/data000-1/data000-1'] ; % image repeated bright forest
DataBlock(9).ImJitterConcat{2} = [RootAnalysisDir,'2017-03-17-0/data000-2/data000-2'] ; % square repeated (no bg - bright forest)
DataBlock(9).ImJitterConcat{3} = [RootAnalysisDir,'2017-03-17-0/data004-5/data004-5'] ; % image and square (bright forest)

% rat
DataBlock(10).BwPath{1} = [RootAnalysisDir,'2017-01-16-0/data006_KR/data006_KR'] ;
DataBlock(10).MovieRepPath{1} = [RootAnalysisDir,'2017-01-16-0/data009-map_KR/data009-map_KR'] ;
DataBlock(10).DsPath{1} = [RootAnalysisDir,'2017-01-16-0/data007-map/data007-map'] ;
DataBlock(10).BwImagePath = [RootImageDir,'2016-12-07-0/4x_oled_mapping_stim_7029.jpg'] ; % BW 15?
DataBlock(10).MoviePath{1} = [RootNatStimMoviesDir,'CatCam/cat_mean117_sd62_0to255.mat'] ; % cat cam
DataBlock(10).MovieStixWidth = 2 ; % movie stixel width
DataBlock(10).MovieFrameInterval = 1 ;% movie frame interval
DataBlock(10).MovieRepFrameNum = 300 ; % movie frame number

% c57/bl6
% may have some trigger issues
DataBlock(11).DsPath{1} = [RootAnalysisDir,'2017-05-04-0/data000/data000'] ;
DataBlock(11).DsPath{2} = [RootAnalysisDir,'2017-05-04-0/data002/data002'] ; % better ds data
DataBlock(11).ImSlide{1} = [RootAnalysisDir,'2017-05-04-0/data001/data001'] ; % forest image, delayed between trigger and stimulus
DataBlock(11).ImSlide{2} = [RootAnalysisDir,'2017-05-04-0/data003/data003'] ; % trees and building image, delayed between trigger and stimulus
DataBlock(11).ImSlide{3} = [RootAnalysisDir,'2017-05-04-0/data004/data004'] ; % grass field image, fixed trig
DataBlock(11).ImSlide{4} = [RootAnalysisDir,'2017-05-04-0/data006/data006'] ; % forest image, fixed trig
DataBlock(11).ImSlideConcat{1} = [RootAnalysisDir,'2017-05-04-0/data002-1/data002-1'] ;
DataBlock(11).ImSlideConcat{2} = [RootAnalysisDir,'2017-05-04-0/data002-3/data002-3'] ; % 
DataBlock(11).ImSlideConcat{3} = [RootAnalysisDir,'2017-05-04-0/data002-4_JC/data002-4_JC'] ; % Best image data
DataBlock(11).ImSlideConcat{4} = [RootAnalysisDir,'2017-05-04-0/data002-6/data002-6'] ; %

% c57/bl6
% not a good recording
% may have some trigger issues
DataBlock(12).DsPath{1} = [RootAnalysisDir,'2017-06-20-0/data000/data000'] ;
DataBlock(12).DsPath{2} = [RootAnalysisDir,'2017-06-20-0/data002/data002'] ; % 
DataBlock(12).ImSlide{1} = [RootAnalysisDir,'2017-06-20-0/data001/data001'] ; % 3 images cycled (building,forest, beach), delayed between trigger and stimulus
DataBlock(12).ImSlideConcat{1} = [RootAnalysisDir,'2017-06-20-0/data000-1/data000-1'] ;
DataBlock(12).ImSlideConcat{2} = [RootAnalysisDir,'2017-06-20-0/data002-1/data002-1'] ;

% c57/bl6
% not a good recording (again)
% may have some trigger issues
DataBlock(13).DsPath{1} = [RootAnalysisDir,'2017-06-23-0/data000/data000'] ;
DataBlock(13).DsPath{2} = [RootAnalysisDir,'2017-06-23-0/data002/data002'] ; % 
DataBlock(13).ImSlide{1} = [RootAnalysisDir,'2017-06-23-0/data001/data001'] ; % 3 images cycled (building,forest, beach), delayed between trigger and stimulus
DataBlock(13).ImSlide{2} = [RootAnalysisDir, '2017-06-23-0/data003/data003'] ; % 13 images - major stim timining variability!
DataBlock(13).ImSlideConcat{1} = [RootAnalysisDir,'2017-06-23-0/data000-1/data000-1'] ;
DataBlock(13).ImSlideConcat{2} = [RootAnalysisDir,'2017-06-23-0/data002-1/data002-1'] ;

% c57/bl6
% great recording 
% may have some trigger issues
DataBlock(14).DsPath{1} = [RootAnalysisDir,'2017-07-21-0/data001/data001'] ; % complete 
DataBlock(14).DsPath{2} = [RootAnalysisDir,'2017-07-21-0/data002/data002'] ; % same directions as the ImageSlide - not uniform across 360 
DataBlock(14).ImSlide{1} = [RootAnalysisDir,'2017-07-21-0/data000/data000'] ; % 3 images cycled (building,forest, beach) in 5 degree increments with 1 at 180 degree
DataBlock(14).ImSlideConcat{1} = [RootAnalysisDir,'2017-07-21-0/data001-0/data001-0'] ;
DataBlock(14).ImSlideDgConcat{1} = [RootAnalysisDir,'2017-07-21-0/data001-0-2/data001-0-2'] ; % 
DataBlock(14).BwPath{1} = [RootAnalysisDir,'2017-07-21-0/data003/data003'] ; % 10-1 

% c57/bl6
% seemed to be good recording
DataBlock(15).DsPath{1} = [RootAnalysisDir,'2017-09-14-0/data002/data002'] ; % probably better ds-id dataset
DataBlock(15).DsPath{2} = [RootAnalysisDir,'2017-09-14-0/data000/data000'] ;
DataBlock(15).DsPath{3} = [RootAnalysisDir,'2017-09-14-0/data005/data005'] ; 
DataBlock(15).DsPath{4} = [RootAnalysisDir,'2017-09-14-0/data003/data003'] ; % 8 dir, 4tp, 4sp, 2 contrasts, 1 reapeat (incomplete)
DataBlock(15).ImSlide{1} = [RootAnalysisDir,'2017-09-14-0/data001/data001'] ; % 6 images, 8 directions, 10 repeats
DataBlock(15).ImSlide{2} = [RootAnalysisDir,'2017-09-14-0/data006/data006'] ; % 6 images, 8 directions, 10 repeats (same as 001 but maybe better rates)
DataBlock(15).ImSlideConcat{1} = [RootAnalysisDir,'2017-09-14-0/data002-1_JC/data002-1_JC'] ;
DataBlock(15).ImSlideConcat{2} = [RootAnalysisDir,'2017-09-14-0/data000-1/data000-1'] ; % 
DataBlock(15).ImSlideConcat{3} = [RootAnalysisDir,'2017-09-14-0/data005-6/data005-6'] ; % 
DataBlock(15).ImSlideDgConcat{1} = [RootAnalysisDir,'2017-09-14-0/data001-0-2/data001-0-2'] ; % 
DataBlock(15).BwPath{1} = [RootAnalysisDir,'2017-09-14-0/data004/data004'] ; % 15-1 

% c57/bl6
% seemed to be ok recording 
% poor DStype pop overlap
DataBlock(16).DsPath{1} = [RootAnalysisDir,'2017-10-12-0/data002/data002'] ; % probably better ds-id dataset
DataBlock(16).DsPath{2} = [RootAnalysisDir,'2017-10-12-0/data000/data000'] ;
DataBlock(16).DsPath{3} = [RootAnalysisDir,'2017-10-12-0/data005/data005'] ; 
DataBlock(16).DsPath{4} = [RootAnalysisDir,'2017-10-12-0/data003/data003'] ; % 8 dir, 4tp, 4sp, 2 contrasts, 1 reapeat (incomplete)
DataBlock(16).ImSlide{1} = [RootAnalysisDir,'2017-10-12-0/data001/data001'] ; % 6 images, 8 directions, 10 repeats
DataBlock(16).ImSlideConcat{1} = [RootAnalysisDir,'2017-10-12-0/data002-1/data002-1'] ;
DataBlock(16).ImJitterConcat{1} = [RootAnalysisDir,'2017-10-12-0/data002-3/data002-3'] ; % one image (1 with some unfinished repeats)
DataBlock(16).ImJitterConcat{2} = [RootAnalysisDir,'2017-10-12-0/data005-6/data005-6'] ; % one image (1 with no repeats - corrected mistake above)
DataBlock(16).BwPath{1} = [RootAnalysisDir,'2017-10-12-0/data004/data004'] ; % 10-2

% c57/bl6
% ok recording - high background but strong response
DataBlock(17).ConcatRevMoviePath{1} = [RootAnalysisDir,'2018-01-26-0/data002-3/data002-3'] ; % cat [movie, movie in reverse]

% c57/bl6
% moving bar data too
% ok recording 
DataBlock(18).DsPath{1} = [RootAnalysisDir,'2018-05-04-0/data000/data000'] ;
DataBlock(18).ImSlide{1} = [RootAnalysisDir,'2018-05-04-0/data001/data001'] ; % 6 images, 8 directions, 10 repeats
DataBlock(18).ImSlideConcat{1} = [RootAnalysisDir,'2018-05-04-0/data000-1_JC/data000-1_JC'] ;

% CBA/C57 wt
% also BW data and longer DG stim at mulitple light levels
DataBlock(19).DsPath{1} = [RootAnalysisDir,'2018-06-01-0/data005/data005'] ;
DataBlock(19).ImJitterConcat{1} = [RootAnalysisDir,'2018-06-01-0/data005-4_JC/data005-4_JC'] ; %image (forest)
DataBlock(19).ImJitterConcat{2} = [RootAnalysisDir,'2018-06-01-0/data005-6/data005-6'] ; %square (no bg - forest)
DataBlock(19).ImJitter{1} = [RootAnalysisDir,'2018-06-01-0/data004/data004'] ; %image (forest)
DataBlock(19).ImJitter{2} = [RootAnalysisDir,'2018-06-01-0/data006/data006'] ; %square 

% CBA/C57 wt
% also stim at mulitple light levels
DataBlock(20).DsPath{1} = [RootAnalysisDir,'2018-07-10-0/data003/data003'] ;
DataBlock(20).DsPath{2} = [RootAnalysisDir,'2018-07-10-0/data007/data007'] ; % very poor ds separation
DataBlock(20).ImJitterConcat{1} = [RootAnalysisDir,'2018-07-10-0/data003-4_JC/data003-4_JC'] ; %image (forest)
DataBlock(20).ImJitterConcat{2} = [RootAnalysisDir,'2018-07-10-0/data003-5/data003-5'] ; %square (no bg - forest)
DataBlock(20).ImJitterConcat{3} = [RootAnalysisDir,'2018-07-10-0/data007-6/data007-6'] ; %image (different forest)
DataBlock(20).ImJitter{1} = [RootAnalysisDir,'2018-07-10-0/data004/data004'] ; %image (forest)
DataBlock(20).ImJitter{2} = [RootAnalysisDir,'2018-07-10-0/data005/data005'] ; %square (no bg - forest)

% c57 wt
% recording also has BW and "direction flows from zucker lab"
DataBlock(21).DsPath{1} = [RootAnalysisDir,'2018-09-26-0/data002/data002'] ;
DataBlock(21).DsPath{2} = [RootAnalysisDir,'2018-09-26-0/data000/data000'] ; % more TF for full speed tuning
DataBlock(21).ImJitterConcat{1} = [RootAnalysisDir,'2018-09-26-0/data002-1/data002-1'] ; %image (forest)
DataBlock(21).ImJitterConcat{2} = [RootAnalysisDir,'2018-09-26-0/data002-3/data002-3'] ; %image (building and trees)
DataBlock(21).ImJitterConcat{3} = [RootAnalysisDir,'2018-09-26-0/data000-1/data000-1'] ; %image (forest) - with different MG
DataBlock(21).ImJitterRepeats = [RootAnalysisDir,'2018-09-26-0/data004/data004'] ; %image (forest)
DataBlock(21).ImJitter{1} = [RootAnalysisDir,'2018-09-26-0/data001/data001'] ; %image (forest)
DataBlock(21).ImJitter{2} = [RootAnalysisDir,'2018-09-26-0/data003/data003'] ; %image (building and trees)



%% Params defaults
Params.empty = [] ;
%ForIgor = [] ;

%% run functions

if perform.MovieRasterPlots ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = MovieRasterPlots(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.MovieSpikeParamFinder ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = MovieSpikeParamFinder(DataBlock, DB, Params) ;
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

if perform.DsAnalysisNatStim ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = DsAnalysisNatStim(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.TformCalculator ;
     for a=1:length(DBset{DBset_id}) ; % for each data block in set
         DB = DBset{DBset_id}(a) ;
         stixSize = 50 ; % assumes a BW 50
         Tform = map_array_from_BW_image(DataBlock(DB).BwImagePath,stixSize,DataBlock(DB).BwMoviePath,DataBlock(DB).ArrayImagePath, DataBlock(DB).BwPath{1}) ;
         save(DataBlock(DB).TformEiPath,'Tform') ;
     end
end
   
if perform.SquareAnalysisV2 ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = SquareAnalysisV2(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.DsAnalysisNatStimV2 ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = DsAnalysisNatStimV2(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.DsAnalysisNatStimV3 ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = DsAnalysisNatStimV3(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.NatStimMappedBw ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = NatStimMappedBw(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end


if perform.ImJitterAnalysis ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = ImJitterAnalysis(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.ImJitterAnalysisV2 ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = ImJitterAnalysisV2(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end


if perform.RepStimNecSufAnalysis ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = RepStimNecSufAnalysis(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.RepStimNecSufAnalysisV2 ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = RepStimNecSufAnalysisV2(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end
   
if perform.ImSliderAnalysis ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = ImSliderAnalysis(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.ImSliderConcatAnalysis ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = ImSliderConcatAnalysis(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.ImSliderConcatAnalysisV2 ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = ImSliderConcatAnalysisV2(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.ImSliderConcatAnalysisV3 ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = ImSliderConcatAnalysisV3(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.ImSliderConcatAnalysisV4 ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = ImSliderConcatAnalysisV4(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.ImSliderConcatAnalysisV5 ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = ImSliderConcatAnalysisV5(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.ImSliderConcatAnalysisV6 ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = ImSliderConcatAnalysisV6(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end


if perform.RepStimClusterAnalysis ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = RepStimClusterAnalysis(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.RevMovConcatAnalyzer ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = RevMovConcatAnalyzer(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

%% population analysis

%% export to hdf5 file

if ForIgor.exportFlag==1 ;
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
%         cd '/Users/jcafaro/Documents/AnalysisFigures/NatStimDs/'  
%         exportStructToHDF5(ForIgor,'DsAnalysisNatStimV3.h5','/')
    else disp('name too long for igor')
    end
end

end % DataBlocks function 


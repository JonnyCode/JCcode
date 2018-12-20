function [DataBlock,Params] = DataBlocks_KO 

% script to analyze amacrine KO experiments (run as script or function)

% 11/11/2014 JC

% functions to perform
perform.KoBwAnalyzer = false ;
perform.KoBwAnalyzer2 = false ;
perform.KoBwAnalyzer3 = false ;
perform.KoBwAnalyzer4 = false ;
perform.KoBwAnalyzer5 = false ;

perform.KoStepAnalyzer = false ;

perform.KoDriftingGratingAnalyzer = false ;

perform.KoSpatialTuningAnalyzer = false ;

perform.KoSpatialTuningAnalyzerNoMapping = false ;
perform.KoSpatialTuningAnalyzerNoMappingV2 = false ;

perform.KoSpatialTuningAnalyzerConcat = false ;
perform.KoSpatialTuningAnalyzerConcatPlusMapping = false ;
perform.KoSpatialTuningAnalyzerConcatPlusMappingV2 = false ;
perform.KoSpatialTuningAnalyzerConcatPlusMappingV3 = false ;
perform.KoSpatialTuningAnalyzerConcatPlusMappingV4 = false ;

perform.KoStepConcatAnalyzer = false ;
perform.KoStepConcatAnalyzerPlusMapping = false ;

perform.DsCellFinder = false ;

perform.KoReverseGratingsAnalyzerConcatPlusMappingV1 = false ;

perform.KoReverseGratingsAnalyzerConcatPlusMappingMultiSetV1 = false ;

perform.KoGratingsAnalyzerConcatPlusMappingMultiSetV1 = false ;

perform.KoBwRepeatConcatAnalyzerPlusMapping = false ;

perform.rfProfileAnalysis = false; 

perform.KoDsConcatAnalysis = false ;

perform.MbDrugAnalysis = false ;

% perform functions on this data block set
DBset_id = 1 ;

% data block sets
DBset{1} = [49] ; % temp analysis
DBset{2} = [3,5,6] ; % VIP-cre mouse with bw stim using map-ei
DBset{3} = [4] ; % VIP-cre mouse with bw stim mapping with Java code only
DBset{4} = [7,8] ; % VIP-cre mouse with ff light step stim using map-ei
DBset{5} = [9,10] ; % VIP-cre mouse with drifiting gratings using map-ei
DBset{6} = [11,12,16:23] ; % Cx57-cre mice +PSEM with DG
DBset{7} = [13,14,15] ; % DG controls
DBset{8} = [21,22,24] ; % Cx57-cre mice +PSEM with BW
DBset{6} = [11,12,16:21] ; % Cx57-cre mice + PSEM with ff steps
DBset{7} = [25] ; % Cx57-icre mice +PSEM with RG
DBset{8} = [26:30] ; % control mice +PSEM with RG
DBset{9} = [31] ; % Cx57-cre mice +low PSEM with RG 
DBset{10} = [32:34] ; % Cx57-cre mice +low PSEM with DG
DBset{11} = [35] ; % DG +/- APV

% Root dircetory for analysis files
%RootAnalysisDir = '/Volumes/Berlioz/Analysis/' ; % Berlioz
RootAnalysisDir = '/Volumes/lab/Experiments/Array/Analysis/' ; % Brahms server
RootMovieDir = '/Volumes/lab/acquisition/movie-xml/' ; % Brahms server xml movie root
RootImageDir = '/Volumes/lab/Experiments/Array/Images/' ; % 

% VIP-cre mouse - bw stim
% Greg thought this one worked
DataBlock(1).PreKo.DataPath = [RootAnalysisDir,'2013-12-22-2/data000-3200-5000s-gdf/data000-3200-5000s'] ; % control
DataBlock(1).PostKo.DataPath = [RootAnalysisDir,'2013-12-22-2/data004/data004'] ; % +psem
DataBlock(1).Wash.DataPath = [RootAnalysisDir,'2013-12-22-2/data006/data006'] ; % wash
DataBlock(1).mapEiFlag = true ; % map by elictrical image

% same retina as db 1 analyzed diffrently 
DataBlock(2).PreKo.DataPath = [RootAnalysisDir,'2013-12-22-2/data000-3200-5000s-gdf/data000-3200-5000s'] ; % control+10
DataBlock(2).PostKo.DataPath = [RootAnalysisDir,'2013-12-22-2/data004-data000-3200-5000s-gdf-map/data004-data000-3200-5000s-gdf-map'] ; % +psem
DataBlock(2).Wash.DataPath = [RootAnalysisDir,'2013-12-22-2/data006/data006'] ; % wash
DataBlock(2).mapEiFlag = false ; % don't map by ei

% same retina as db 1 but organized differently for KoBwAnalyzer2
DataBlock(3).DataPath{1} = [RootAnalysisDir,'2013-12-22-2/data000-3200-5000s-gdf/data000-3200-5000s'] ; % control
DataBlock(3).DataPath{1} = [RootAnalysisDir,'2013-12-22-2/data000-3200-5000s-JC/data000-3200-5000s'] ; % control with "post analysis" group
DataBlock(3).DataPath{2} = [RootAnalysisDir,'2013-12-22-2/data004/data004'] ; % +PSEM
DataBlock(3).DataPath{3} = [RootAnalysisDir,'2013-12-22-2/data006/data006'] ; % wash
DataBlock(3).mapEiFlag = true ;

% same data analysis as db 2 but organized for KoBwAnalyzer2
DataBlock(4).DataPath{1} = [RootAnalysisDir,'2013-12-22-2/data000-3200-5000s-gdf/data000-3200-5000s'] ; % control
DataBlock(4).DataPath{2} = [RootAnalysisDir,'2013-12-22-2/data004-data000-3200-5000s-gdf-map/data004-data000-3200-5000s-gdf-map'] ; % +PSEM
DataBlock(4).DataPath{3} = [RootAnalysisDir,'2013-12-22-2/data006/data006'] ; % wash
DataBlock(4).mapEiFlag = false ;

% VIP cre mouse - bw stim
DataBlock(5).DataPath{1} = [RootAnalysisDir,'2014-12-02-0/data002/data002'] ; % control
DataBlock(5).DataPath{1} = [RootAnalysisDir,'2014-12-02-0/data002-JC/data002'] ; % control with "post analysis" group
DataBlock(5).DataPath{2} = [RootAnalysisDir,'2014-12-02-0/data004/data004'] ; % +PSEM
DataBlock(5).DataPath{3} = [RootAnalysisDir,'2014-12-02-0/data008/data008'] ; % wash
DataBlock(5).mapEiFlag = true ;

% VIP cre mouse - bw stim
% this data is missing first trigger - cant get STAs with vision-auto-sta !!!!!!!!
DataBlock(6).DataPath{1} = [RootAnalysisDir,'2014-12-06-0/data001/data001'] ; % control
DataBlock(6).DataPath{2} = [RootAnalysisDir,'2014-12-06-0/data004/data004'] ; % +PSAM (psem?)
DataBlock(6).DataPath{3} = [RootAnalysisDir,'2014-12-06-0/data007/data007'] ; % wash
DataBlock(6).mapEiFlag = true ;

% VIP cre mouse - full field light steps gwgb
% same retina as db 5 bw data
DataBlock(7).DataPath{1} = [RootAnalysisDir,'2014-12-02-0/data002/data002'] ; % control white noise
DataBlock(7).DataPath{2} = [RootAnalysisDir,'2014-12-02-0/data000/data000'] ; % control 
DataBlock(7).DataPath{3} = [RootAnalysisDir,'2014-12-02-0/data003/data003'] ; % +PSEM early
DataBlock(7).DataPath{4} = [RootAnalysisDir,'2014-12-02-0/data006/data006'] ; % +PSEM late
DataBlock(7).DataPath{5} = [RootAnalysisDir,'2014-12-02-0/data007/data007'] ; % wash early
DataBlock(7).mapEiFlag = true ;

% VIP cre mouse - full field light steps gwgb
% same retina as db 6 bw data
% notes mention a change in light step during psem perfusion and wash - compare early and late
DataBlock(8).DataPath{1} = [RootAnalysisDir,'2014-12-06-0/data001/data001'] ; % control white noise
DataBlock(8).DataPath{2} = [RootAnalysisDir,'2014-12-06-0/data002/data002'] ; % control  
DataBlock(8).DataPath{3} = [RootAnalysisDir,'2014-12-06-0/data003/data003'] ; % + PSEM early and late 
DataBlock(8).DataPath{4} = [RootAnalysisDir,'2014-12-06-0/data006/data006'] ; % wash

% VIP cre mouse - drifting gratings
% same retina as db 5
DataBlock(9).DataPath{1} = [RootAnalysisDir,'2014-12-02-0/data002/data002'] ; % control white noise
DataBlock(9).DataPath{2} = [RootAnalysisDir,'2014-12-02-0/data001/data001'] ; % control
DataBlock(9).DataPath{3} = [RootAnalysisDir,'2014-12-02-0/data005/data005'] ; % +PSEM
DataBlock(9).DataPath{4} = [RootAnalysisDir,'2014-12-02-0/data009/data009'] ; % wash
DataBlock(9).StimPath{2} = [RootAnalysisDir,'2014-12-02-0/stimuli/s01'] ; % control
DataBlock(9).StimPath{3} = [RootAnalysisDir,'2014-12-02-0/stimuli/s05'] ; % +PSEM
DataBlock(9).StimPath{4} = [RootAnalysisDir,'2014-12-02-0/stimuli/s09'] ; % wash
DataBlock(9).mapEiFlag = true ;

% VIP cre mouse - drifting gratings
% same retina as db 6
DataBlock(10).DataPath{1} = [RootAnalysisDir,'2014-12-06-0/data001/data001'] ; % control white noise
DataBlock(10).DataPath{2} = [RootAnalysisDir,'2014-12-06-0/data000/data000'] ; % control
DataBlock(10).DataPath{3} = [RootAnalysisDir,'2014-12-06-0/data005/data005'] ; % +PSAM (psem?)
DataBlock(10).DataPath{4} = [RootAnalysisDir,'2014-12-06-0/data008/data008'] ; % wash
DataBlock(10).StimPath{2} = [RootAnalysisDir,'2014-12-06-0/stimuli/s00'] ; % control
DataBlock(10).StimPath{3} = [RootAnalysisDir,'2014-12-06-0/stimuli/s05'] ; % +PSAM (psem?)
DataBlock(10).StimPath{4} = [RootAnalysisDir,'2014-12-06-0/stimuli/s08'] ; % wash
DataBlock(10).mapEiFlag = true ;

% Cx57 cre mouse - drifiting gratings
DataBlock(11).BwPath{1} = [RootAnalysisDir,'2015-11-12-0/data004/data004'] ; % control white noise (BW 20-1)
DataBlock(11).DgPath{1} = [RootAnalysisDir,'2015-11-12-0/data003/data003'] ; % control 25% contrast
DataBlock(11).DgPath{2} = [RootAnalysisDir,'2015-11-12-0/data007/data007'] ; % +PSEM 25% conrast
DataBlock(11).DgPath{3} = [RootAnalysisDir,'2015-11-12-0/data011/data011'] ; % wash 25% contrast
DataBlock(11).DgPath{4} = [RootAnalysisDir,'2015-11-12-0/data012/data012'] ; % wash 100% contrast
DataBlock(11).DgConcat = [RootAnalysisDir,'2015-11-12-0/data003-7-11/data003-7-11'] ; % wash 100% contrast
DataBlock(11).FfPulse{1} = [RootAnalysisDir,'2015-11-12-0/data005/data005'] ; % +PSEM ~400s
DataBlock(11).FfPulse{2} = [RootAnalysisDir,'2015-11-12-0/data009/data009'] ; % -PSEM ~400s
DataBlock(11).FfPulseConcat = [RootAnalysisDir,'2015-11-12-0/data005-9/data005-9'] ; % -PSEM ~400s
DataBlock(11).FfPulseConcatDrugTimes = [400, 1900] ;
DataBlock(11).mapEiFlag = true ;
DataBlock(11).DataPath{1} = [RootAnalysisDir,'2015-11-12-0/data004/data004'] ; % control white noise (BW 20-1)
DataBlock(11).DataPath{2} = [RootAnalysisDir,'2015-11-12-0/data006/data006'] ; % +PSEM (BW 20-1)
DataBlock(11).DataPath{3} = [RootAnalysisDir,'2015-11-12-0/data010/data010'] ; % -PSEM (BW 20-1)
% DataBlock(11).DataPath{1} = [RootAnalysisDir,'2015-11-12-0/data001/data001'] ; % control white noise (BW 10-4)
% DataBlock(11).DataPath{2} = [RootAnalysisDir,'2015-11-12-0/data008/data008'] ; % +PSEM (BW 10-4)



% Cx57 cre mouse - drifiting gratings
% 50/50 mirror with 0.3 NDF in front (attenuation ~2)
%DataBlock(12).BwPath{1} = [RootAnalysisDir,'2015-11-17-0/data006/data006'] ; % control white noise (BW 15-2)
DataBlock(12).BwPath{1} = [RootAnalysisDir,'2015-11-17-0/data006_RS/data006_RS'] ; % control white noise (BW 15-2)
DataBlock(12).DgPath{1} = [RootAnalysisDir,'2015-11-17-0/data000/data000'] ; % control 12,24,48% contrast
DataBlock(12).DgPath{2} = [RootAnalysisDir,'2015-11-17-0/data002/data002'] ; % +PSEM 12,24,48% contrast
DataBlock(12).DgPath{3} = [RootAnalysisDir,'2015-11-17-0/data005/data005'] ; % wash 12,24,48% contrast
DataBlock(12).DgConcat = [RootAnalysisDir,'2015-11-17-0/data000-2-5_RS/data000-2-5_RS'] ; % wash 12,24,48% contrast
DataBlock(12).FfPulse{1} = [RootAnalysisDir,'2015-11-17-0/data001/data001'] ; % +PSEM at ~300s 48% contrast
DataBlock(12).FfPulse{2} = [RootAnalysisDir,'2015-11-17-0/data003-4/data003-4'] ; % -PSEM at ~228s 48% contrast
DataBlock(12).FfPulseConcat = [RootAnalysisDir,'2015-11-17-0/data001-3-4_RS/data001-3-4_RS'] ; % -PSEM at ~228s 48% contrast
DataBlock(12).FfPulseConcatDrugTimes = [300, 1700] ; % [+psem,-psem]NOT EXACT YET!!!!
DataBlock(12).DsPath = [RootAnalysisDir,'2015-11-17-0/data007/data007'] ; % moving bars in multidirections for finding DS cells
DataBlock(12).mapEiFlag = true ;

% CBA mouse - drifiting gratings - control for Cx-PSAM experiments (no PSEM in wt CBA)
% 50/50 mirror with 0.3 NDF in front (attenuation ~2)
DataBlock(13).BwPath{1} = [RootAnalysisDir,'2016-01-25-0/data005/data005'] ; % control white noise (BW 15-2)
DataBlock(13).DgPath{1} = [RootAnalysisDir,'2016-01-25-0/data000/data000'] ; % control 12,24,48% contrast
DataBlock(13).DgPath{2} = [RootAnalysisDir,'2016-01-25-0/data002/data002'] ; % different Ames solution 12,24,48% contrast
DataBlock(13).DgPath{3} = [RootAnalysisDir,'2016-01-25-0/data004/data004'] ; % wash 12,24,48% contrast
DataBlock(13).DgConcat = [RootAnalysisDir,'2016-01-25-0/data000-2-4/data000-2-4'] ; % wash 12,24,48% contrast
DataBlock(13).FfPulse{1} = [RootAnalysisDir,'2016-01-25-0/data001/data001'] ; % +new Ames at ~320s 48% contrast
DataBlock(13).FfPulse{2} = [RootAnalysisDir,'2016-01-25-0/data003/data003'] ; % -new Ames at ~320s 48% contrast
DataBlock(13).FfPulseConcat = [RootAnalysisDir,'2016-01-25-0/data001-3/data001-3'] ; % 
DataBlock(13).FfPulseConcatDrugTimes = [320, 1720] ; % 
DataBlock(13).DsPath = [RootAnalysisDir,'2016-01-25-0/data006/data006'] ; % moving bars in multidirections for finding DS cells
DataBlock(13).mapEiFlag = true ;

% CBA mouse - drifiting gratings - control for Cx-PSAM experiments (+/-PSEM wt CBA)
% 50/50 mirror with 0.3 NDF in front (attenuation ~2)
DataBlock(14).BwPath{1} = [RootAnalysisDir,'2016-02-03-0/data005/data005'] ; % control white noise (BW 15-2)
DataBlock(14).DgPath{1} = [RootAnalysisDir,'2016-02-03-0/data000/data000'] ; % control 12,24,48% contrast
DataBlock(14).DgPath{2} = [RootAnalysisDir,'2016-02-03-0/data002/data002'] ; % +PSEM solution 12,24,48% contrast
DataBlock(14).DgPath{3} = [RootAnalysisDir,'2016-02-03-0/data004/data004'] ; % wash 12,24,48% contrast
DataBlock(14).DgConcat = [RootAnalysisDir,'2016-02-03-0/data000-2-4/data000-2-4'] ; % wash 12,24,48% contrast
DataBlock(14).FfPulse{1} = [RootAnalysisDir,'2016-02-03-0/data001/data001'] ; % +PSEM at ~s 48% contrast
DataBlock(14).FfPulse{2} = [RootAnalysisDir,'2016-02-03-0/data003/data003'] ; % -PSEM at ~s 48% contrast
DataBlock(14).FfPulseConcat = [RootAnalysisDir,'2016-02-03-0/data001-3/data001-3'] ; % 
DataBlock(14).FfPulseConcatDrugTimes = [300, 1725] ; % 
DataBlock(14).DsPath = [RootAnalysisDir,'2016-02-03-0/data006/data006'] ; % moving bars in multidirections for finding DS cells
DataBlock(14).mapEiFlag = true ;

% CBA mouse - drifiting gratings - control for Cx-PSAM experiments (+/-PSEM wt CBA)
% 50/50 mirror with 0.3 NDF in front (attenuation ~2)
% THIS ENTIRE DATASET IS MISSING ALL TRIGGERS!!!! 
DataBlock(15).BwPath{1} = [RootAnalysisDir,'2016-02-05-0/data005/data005'] ; % control white noise (BW 15-2) 
DataBlock(15).DgPath{1} = [RootAnalysisDir,'2016-02-05-0/data000/data000'] ; % control 12,24,48% contrast
DataBlock(15).DgPath{2} = [RootAnalysisDir,'2016-02-05-0/data002/data002'] ; % +PSEM solution 12,24,48% contrast
DataBlock(15).DgPath{3} = [RootAnalysisDir,'2016-02-05-0/data004/data004'] ; % wash 12,24,48% contrast
DataBlock(15).DgConcat = [RootAnalysisDir,'2016-02-05-0/data000-2-4/data000-2-4'] ; % wash 12,24,48% contrast
DataBlock(15).FfPulse{1} = [RootAnalysisDir,'2016-02-05-0/data001/data001'] ; % +PSEM at ~s 48% contrast
DataBlock(15).FfPulse{2} = [RootAnalysisDir,'2016-02-05-0/data003/data003'] ; % -PSEM at >s 48% contrast
DataBlock(15).FfPulseConcat = [RootAnalysisDir,'2016-02-05-0/data001-3/data001-3'] ; % 
DataBlock(15).FfPulseConcatDrugTimes = [] ; % 
DataBlock(15).DsPath = [RootAnalysisDir,'2016-02-05-0/data006/data006'] ; % moving bars in multidirections for finding DS cells
DataBlock(15).mapEiFlag = true ;

% Cx57-Cre + AAV7m8-hEF1a-flex-PSAM - drifiting gratings  (+/-PSEM)
% 50/50 mirror with 0.3 NDF in front (attenuation ~2)
% transfection looked good, flow rate may have been low (see notes)
DataBlock(16).BwPath{1} = [RootAnalysisDir,'2016-03-09-0/data005/data005'] ; % control white noise (BW 15-2) 
DataBlock(16).DgPath{1} = [RootAnalysisDir,'2016-03-09-0/data000/data000'] ; % control 12,24,48% contrast
DataBlock(16).DgPath{2} = [RootAnalysisDir,'2016-03-09-0/data002/data002'] ; % +PSEM solution 12,24,48% contrast
DataBlock(16).DgPath{3} = [RootAnalysisDir,'2016-03-09-0/data004/data004'] ; % wash 12,24,48% contrast
DataBlock(16).DgConcat = [RootAnalysisDir,'2016-03-09-0/data000-2-4/data000-2-4'] ; % wash 12,24,48% contrast
DataBlock(16).FfPulse{1} = [RootAnalysisDir,'2016-03-09-0/data001/data001'] ; % +PSEM at ~340s 48% contrast
DataBlock(16).FfPulse{2} = [RootAnalysisDir,'2016-03-09-0/data003/data003'] ; % -PSEM at >320s 48% contrast
DataBlock(16).FfPulseConcat = [RootAnalysisDir,'2016-03-09-0/data001-3/data001-3'] ; % 
DataBlock(16).FfPulseConcatDrugTimes = [340, 1720] ; % 
DataBlock(16).DsPath = [RootAnalysisDir,'2016-03-09-0/data006/data006'] ; % moving bars in multidirections for finding DS cells
DataBlock(16).mapEiFlag = true ;

% Cx57-Cre + AAV7m8-hEF1a-flex-PSAM - drifiting gratings  (+/-PSEM)
% 50/50 mirror with 0.3 NDF in front (attenuation ~2)
% transfection looked good
% changed dg contrast from previous experiments
DataBlock(17).BwPath{1} = [RootAnalysisDir,'2016-03-10-0/data005/data005'] ; % control white noise (BW 15-2) 
DataBlock(17).DgPath{1} = [RootAnalysisDir,'2016-03-10-0/data000/data000'] ; % control 24,48,52% contrast
DataBlock(17).DgPath{2} = [RootAnalysisDir,'2016-03-10-0/data002/data002'] ; % +PSEM solution 24,48,52% contrast
DataBlock(17).DgPath{3} = [RootAnalysisDir,'2016-03-10-0/data004/data004'] ; % wash 24,48,52% contrast
DataBlock(17).DgConcat = [RootAnalysisDir,'2016-03-10-0/data000-2-4/data000-2-4'] ; % wash 12,24,48% contrast
DataBlock(17).FfPulse{1} = [RootAnalysisDir,'2016-03-10-0/data001/data001'] ; % +PSEM at ~340s 48% contrast
DataBlock(17).FfPulse{2} = [RootAnalysisDir,'2016-03-10-0/data003/data003'] ; % -PSEM at >320s 48% contrast
DataBlock(17).FfPulseConcat = [RootAnalysisDir,'2016-10-10-0/data001-3/data001-3'] ; % 
DataBlock(17).FfPulseConcatDrugTimes = [340, 1720] ; % 
DataBlock(17).DsPath = [RootAnalysisDir,'2016-03-10-0/data006/data006'] ; % moving bars in multidirections for finding DS cells
DataBlock(17).mapEiFlag = true ;

% Cx57-Cre + AAV7m8-hEF1a-flex-PSAM - drifiting gratings  (+/-PSEM)
% 50/50 mirror with 0.3 NDF in front (attenuation ~2)
% transfection looked very good
% used hgiher dg contrast as in 2016-3-10 experiments
DataBlock(18).BwPath{1} = [RootAnalysisDir,'2016-03-14-0/data005/data005'] ; % control white noise (BW 15-2) 
DataBlock(18).DgPath{1} = [RootAnalysisDir,'2016-03-14-0/data000/data000'] ; % control 24,48,52% contrast
DataBlock(18).DgPath{2} = [RootAnalysisDir,'2016-03-14-0/data002/data002'] ; % +PSEM solution 24,48,52% contrast
DataBlock(18).DgPath{3} = [RootAnalysisDir,'2016-03-14-0/data004/data004'] ; % wash 24,48,52% contrast
DataBlock(18).DgConcat = [RootAnalysisDir,'2016-03-14-0/data000-2-4/data000-2-4'] ; % wash 12,24,48% contrast
DataBlock(18).FfPulse{1} = [RootAnalysisDir,'2016-03-14-0/data001/data001'] ; % +PSEM at ~320s 48% contrast
DataBlock(18).FfPulse{2} = [RootAnalysisDir,'2016-03-14-0/data003/data003'] ; % -PSEM at >320s 48% contrast
DataBlock(18).FfPulseConcat = [RootAnalysisDir,'2016-03-14-0/data001-3/data001-3'] ; % 
DataBlock(18).FfPulseConcatDrugTimes = [320, 1720] ; % 
DataBlock(18).DsPath = [RootAnalysisDir,'2016-03-14-0/data006/data006'] ; % moving bars in multidirections for finding DS cells
DataBlock(18).mapEiFlag = true ;

% Cx57-Cre + AAV7m8-hEF1a-flex-PSAM - drifiting gratings  (+/-PSEM)
% 50/50 mirror with 0.3 NDF in front (attenuation ~2)
% no gfp+ cells were visible - transfection may have failed
% used lower dg contrast as in 2016-11-17 experiments
DataBlock(19).BwPath{1} = [RootAnalysisDir,'2016-03-16-0/data000/data000'] ; % control white noise (BW 15-2) 
DataBlock(19).DgPath{1} = [RootAnalysisDir,'2016-03-16-0/data002/data002'] ; % control 12,24,48% contrast
DataBlock(19).DgPath{2} = [RootAnalysisDir,'2016-03-16-0/data004/data004'] ; % +PSEM solution 12,24,48% contrast
DataBlock(19).DgPath{3} = [RootAnalysisDir,'2016-03-16-0/data006/data006'] ; % wash 12,24,48% contrast
DataBlock(19).DgConcat = [RootAnalysisDir,'2016-03-16-0/data002-4-6/data002-4-6'] ; % 12,24,48% contrast
DataBlock(19).FfPulse{1} = [RootAnalysisDir,'2016-03-16-0/data003/data003'] ; % +PSEM at ~320s 48% contrast
DataBlock(19).FfPulse{2} = [RootAnalysisDir,'2016-03-16-0/data005/data005'] ; % -PSEM at >320s 48% contrast
DataBlock(19).FfPulseConcat = [RootAnalysisDir,'2016-03-16-0/data003-5/data003-5'] ; % 
DataBlock(19).FfPulseConcatDrugTimes = [320, 1720] ; % 
DataBlock(19).DsPath = [RootAnalysisDir,'2016-03-16-0/data001/data001'] ; % moving bars in multidirections for finding DS cells
DataBlock(19).mapEiFlag = true ;

% Cx57-Cre + AAV7m8-CMV-flex-PSAM - drifiting gratings  (+/-PSEM)
% 50/50 mirror with 0.3 NDF in front (attenuation ~2)
% used lower dg contrast as in 2016-11-17 experiments
% transfection looked ok
DataBlock(20).BwPath{1} = [RootAnalysisDir,'2016-03-18-0/data000/data000'] ; % control white noise (BW 15-2) 
DataBlock(20).DgPath{1} = [RootAnalysisDir,'2016-03-18-0/data002/data002'] ; % control 12,24,48% contrast
DataBlock(20).DgPath{2} = [RootAnalysisDir,'2016-03-18-0/data005/data005'] ; % +PSEM solution 12,24,48% contrast
DataBlock(20).DgPath{3} = [RootAnalysisDir,'2016-03-18-0/data007/data007'] ; % wash 12,24,48% contrast
DataBlock(20).DgConcat = [RootAnalysisDir,'2016-03-18-0/data002-5-7/data002-5-7'] ; % 12,24,48% contrast
DataBlock(20).FfPulse{1} = [RootAnalysisDir,'2016-03-18-0/data003/data003'] ; % +PSEM at ~330s 48% contrast
DataBlock(20).FfPulse{2} = [RootAnalysisDir,'2016-03-18-0/data006/data006'] ; % -PSEM at >215s 48% contrast
DataBlock(20).FfPulseConcat = [RootAnalysisDir,'2016-03-18-0/data003-6/data003-6'] ; % 
DataBlock(20).FfPulseConcatDrugTimes = [330, 1615] ; % 
DataBlock(20).DsPath = [RootAnalysisDir,'2016-03-18-0/data001/data001'] ; % moving bars in multidirections for finding DS cells
DataBlock(20).mapEiFlag = true ;

% Cx57-Cre + AAV7m8-hEF1a-flex-PSAM - drifiting gratings  (+/-PSEM)
% 50/50 mirror with 0.3 NDF in front (attenuation ~2)
% used lower dg contrast as in 2016-11-17 experiments
% no visible gfp+ cells
DataBlock(21).BwPath{1} = [RootAnalysisDir,'2016-03-21-0/data000/data000'] ; % control white noise (BW 15-2) 
DataBlock(21).DgPath{1} = [RootAnalysisDir,'2016-03-21-0/data002/data002'] ; % control 12,24,48% contrast
DataBlock(21).DgPath{2} = [RootAnalysisDir,'2016-03-21-0/data004/data004'] ; % +PSEM solution 12,24,48% contrast
DataBlock(21).DgPath{3} = [RootAnalysisDir,'2016-03-21-0/data006/data006'] ; % wash 12,24,48% contrast
DataBlock(21).DgConcat = [RootAnalysisDir,'2016-03-21-0/data002-4-6/data002-4-6'] ; % 12,24,48% contrast
DataBlock(21).FfPulse{1} = [RootAnalysisDir,'2016-03-21-0/data003/data003'] ; % +PSEM at ~330s 48% contrast
DataBlock(21).FfPulse{2} = [RootAnalysisDir,'2016-03-21-0/data006/data006'] ; % -PSEM at >215s 48% contrast
DataBlock(21).FfPulseConcat = [RootAnalysisDir,'2016-03-21-0/data003-6/data003-6'] ; % 
DataBlock(21).FfPulseConcatDrugTimes = [330, 1615] ; % 
DataBlock(21).DsPath = [RootAnalysisDir,'2016-03-21-0/data001/data001'] ; % moving bars in multidirections for finding DS cells
DataBlock(21).DataPath{1} = [RootAnalysisDir,'2016-03-21-0/data007-0-1800/data007-0-1800'] ; % BW 30-2 control (post previous PSEM-wash)
DataBlock(21).DataPath{2} = [RootAnalysisDir,'2016-03-21-0/data007-2400-4000/data007-2400-4000'] ; % BW 30-2 +PSEM
DataBlock(21).DataPath{3} = [RootAnalysisDir,'2016-03-21-0/data007-4600-6800/data007-4600-6800'] ; %BW 30-2 wash
DataBlock(21).BwMoviePath = [RootMovieDir,'BW-30-2-0.48-11111-20x20-60.35.xml'] ;
DataBlock(21).mapEiFlag = true ;

% Cx57-Cre + AAV7m8-hEF1a-flex-PSAM 2016-03-23 
% Did not get sufficient array coverage to make recording worthwhile
% transfection looked good

% Cx57-Cre + AAV7m8-hEF1a-flex-PSAM - BW  (+/-PSEM)
% 50/50 mirror with 0.3 NDF in front (attenuation ~2)
% transfection looked good
DataBlock(22).BwPath{1} = [RootAnalysisDir,'2016-03-25-0/data000/data000'] ; % control white noise (BW 15-2) 
DataBlock(22).DataPath{1} = [RootAnalysisDir,'2016-03-25-0/data001-900-3000/data001-900-3000'] ; % BW 30-2control (post previous PSEM-wash)
DataBlock(22).DataPath{2} = [RootAnalysisDir,'2016-03-25-0/data001-3900-6000/data001-3900-6000'] ; % BW 30-2 +PSEM
DataBlock(22).DataPath{3} = [RootAnalysisDir,'2016-03-25-0/data001-6900-9000/data001-6900-9000'] ; % BW 30-2 wash
DataBlock(22).DsPath = [RootAnalysisDir,'2016-03-25-0/data002/data002'] ; % moving bars in multidirections for finding DS cells
DataBlock(22).mapEiFlag = true ;
DataBlock(22).BwMoviePath = [RootMovieDir,'BW-30-2-0.48-11111-20x20-60.35.xml'] ;

% Cx57-Cre + AAV7m8-hEF1a-flex-PSAM 2016-03-30 
% Electrical noise prevented recording
% transfection looked good

% Cx57-Cre + AAV7m8-hEF1a-flex-PSAM 2016-04-04 
% Did not get clean mount 
% many gfp+ cells were seen but transfection looked weak

% Cx57-Cre + AAV7m8-2A-hEF1a-flex-PSAM - drifiting gratings  (+/-PSEM)
% 50/50 mirror with 0.3 NDF in front (attenuation ~2)
% td-tom made it difficult to assess transfection 
% used new perfusion switch system for the first time 
DataBlock(23).BwPath{1} = [RootAnalysisDir,'2016-04-20-0/data000/data000'] ; % control white noise (BW 15-2) 
DataBlock(23).DgPath{1} = [RootAnalysisDir,'2016-04-20-0/data002/data002'] ; % control 12,24,48% contrast
DataBlock(23).DgPath{2} = [RootAnalysisDir,'2016-04-20-0/data004/data004'] ; % +PSEM solution 12,24,48% contrast
DataBlock(23).DgPath{3} = [RootAnalysisDir,'2016-04-20-0/data006/data006'] ; % wash 12,24,48% contrast
DataBlock(23).DgConcat = [RootAnalysisDir,'2016-04-20-0/data002-4-6/data002-4-6'] ; % 12,24,48% contrast
DataBlock(23).FfPulse{1} = [RootAnalysisDir,'2016-04-20-0/data003/data003'] ; % +PSEM at ~340s 48% contrast
DataBlock(23).FfPulse{2} = [RootAnalysisDir,'2016-04-20-0/data005/data005'] ; % -PSEM at >340s 48% contrast
DataBlock(23).FfPulseConcat = [RootAnalysisDir,'2016-04-20-0/data003-5/data003-5'] ; % 
DataBlock(23).FfPulseConcatDrugTimes = [340, 1740] ; % 
DataBlock(23).DsPath = [RootAnalysisDir,'2016-04-20-0/data001/data001'] ; % moving bars in multidirections for finding DS cells
DataBlock(23).mapEiFlag = true ;

% Cx57-Cre + AAV7m8-hEF1a-flex-PSAM - BW  (+/-PSEM)
% 50/50 mirror with 0.3 NDF in front (attenuation ~2)
% td-tom made it difficult to assess transfection 
DataBlock(24).BwPath{1} = [RootAnalysisDir,'2016-04-21-0/data000/data000'] ; % control white noise (BW 15-2) 
DataBlock(24).DataPath{1} = [RootAnalysisDir,'2016-04-21-0/data001-900-3000/data001-900-3000'] ; % BW 30-2control (post previous PSEM-wash)
DataBlock(24).DataPath{2} = [RootAnalysisDir,'2016-04-21-0/data001-3900-6000/data001-3900-6000'] ; % BW 30-2 +PSEM
DataBlock(24).DataPath{3} = [RootAnalysisDir,'2016-04-21-0/data001-6900-9000/data001-6900-9000'] ; % BW 30-2 wash
DataBlock(24).DsPath = [RootAnalysisDir,'2016-04-21-0/data002/data002'] ; % moving bars in multidirections for finding DS cells
DataBlock(24).mapEiFlag = true ;
DataBlock(24).BwMoviePath = [RootMovieDir,'BW-30-2-0.48-11111-20x20-60.35.xml'] ;

% Cx57-Cre + AAV7m8-2A-hEF1a-flex-PSAM - contrast reversing gratings  (+/-PSEM)
% 50/50 mirror with 0.3 NDF in front (attenuation ~2)
% td-tom made it difficult to assess transfection 
DataBlock(25).BwPath{1} = [RootAnalysisDir,'2016-05-04-0/data000/data000'] ; % control white noise (BW 15-2) 
DataBlock(25).RgPath{1} = [RootAnalysisDir,'2016-05-04-0/data002/data002'] ; % control 
DataBlock(25).RgPath{2} = [RootAnalysisDir,'2016-05-04-0/data004/data004'] ; % +PSEM solution 
DataBlock(25).RgPath{3} = [RootAnalysisDir,'2016-05-04-0/data006/data006'] ; % wash 
DataBlock(25).RgConcat = [RootAnalysisDir,'2016-05-04-0/data002-4-6/data002-4-6'] ; % 
DataBlock(25).FfPulse{1} = [RootAnalysisDir,'2016-05-04-0/data003/data003'] ; % +PSEM at ~730
DataBlock(25).FfPulse{2} = [RootAnalysisDir,'2016-05-04-0/data005/data005'] ; % -PSEM at ~350 s
DataBlock(25).FfPulseConcat = [RootAnalysisDir,'2016-05-04-0/data003-5/data003-5'] ; % 
DataBlock(25).FfPulseConcatDrugTimes = [730,1742] ; % 
DataBlock(25).DsPath = [RootAnalysisDir,'2016-05-04-0/data001/data001'] ; % moving bars in multidirections for finding DS cells
DataBlock(25).mapEiFlag = true ;

% Control animal, CAUTION: may be Wt, Cx57 het, or Cx57 homo (which is Cx57 KO)
% 50/50 mirror with 0.3 NDF in front (attenuation ~2)
% small spikes, run down - not a great recording
DataBlock(26).BwPath{1} = [RootAnalysisDir,'2016-05-27-0/data000/data000'] ; % control white noise (BW 15-2) 
DataBlock(26).RgPath{1} = [RootAnalysisDir,'2016-05-27-0/data001/data001'] ; % control 
DataBlock(26).RgPath{2} = [RootAnalysisDir,'2016-05-27-0/data003/data003'] ; % +PSEM solution 
DataBlock(26).RgPath{3} = [RootAnalysisDir,'2016-05-27-0/data005/data005'] ; % wash 
DataBlock(26).RgConcat = [RootAnalysisDir,'2016-05-27-0/data001-3-5/data001-3-5'] ; % 
DataBlock(26).FfPulse{1} = [RootAnalysisDir,'2016-05-27-0/data002/data002'] ; % +PSEM at ~380
DataBlock(26).FfPulse{2} = [RootAnalysisDir,'2016-05-27-0/data004/data004'] ; % -PSEM at ~350 s
DataBlock(26).FfPulseConcat = [RootAnalysisDir,'2016-05-27-0/data002-4/data002-4'] ; % 
DataBlock(26).FfPulseConcatDrugTimes = [380,1750] ; % 
DataBlock(26).mapEiFlag = true ;

% Control animal, CAUTION: may be Wt, Cx57 het, or Cx57 homo (which is Cx57 KO)
% 50/50 mirror with 0.3 NDF in front (attenuation ~2)
% seemed to have run-down
DataBlock(27).BwPath{1} = [RootAnalysisDir,'2016-06-02-0/data000/data000'] ; % control white noise (BW 15-2) 
DataBlock(27).RgPath{1} = [RootAnalysisDir,'2016-06-02-0/data002/data002'] ; % control 
DataBlock(27).RgPath{2} = [RootAnalysisDir,'2016-06-02-0/data004/data004'] ; % +PSEM solution 
DataBlock(27).RgPath{3} = [RootAnalysisDir,'2016-06-02-0/data006/data006'] ; % wash 
DataBlock(27).RgConcat = [RootAnalysisDir,'2016-06-02-0/data002-4-6/data002-4-6'] ; % 
DataBlock(27).FfPulse{1} = [RootAnalysisDir,'2016-06-02-0/data003/data003'] ; % +PSEM at ~330
DataBlock(27).FfPulse{2} = [RootAnalysisDir,'2016-06-02-0/data005/data005'] ; % -PSEM at ~340 s
DataBlock(27).FfPulseConcat = [RootAnalysisDir,'2016-06-02-0/data003-5/data003-5'] ; % 
DataBlock(27).FfPulseConcatDrugTimes = [330,1740] ; % 
DataBlock(27).DsPath = [RootAnalysisDir,'2016-06-02-0/data001/data001'] ; % moving bars in multidirections for finding DS cells
DataBlock(27).mapEiFlag = true ;

% Control animal, CAUTION: may be Wt, Cx57 het, or Cx57 homo (which is Cx57 KO)
% 50/50 mirror with 0.3 NDF in front (attenuation ~2)
% seemed to stable
DataBlock(28).BwPath{1} = [RootAnalysisDir,'2016-06-08-0/data000/data000'] ; % control white noise (BW 15-2) 
DataBlock(28).RgPath{1} = [RootAnalysisDir,'2016-06-08-0/data002/data002'] ; % control 
DataBlock(28).RgPath{2} = [RootAnalysisDir,'2016-06-08-0/data004/data004'] ; % +PSEM solution 
DataBlock(28).RgPath{3} = [RootAnalysisDir,'2016-06-08-0/data006/data006'] ; % wash 
DataBlock(28).RgConcat = [RootAnalysisDir,'2016-06-08-0/data002-4-6/data002-4-6'] ; % 
DataBlock(28).FfPulse{1} = [RootAnalysisDir,'2016-06-08-0/data003/data003'] ; % +PSEM at ~360
DataBlock(28).FfPulse{2} = [RootAnalysisDir,'2016-06-08-0/data005/data005'] ; % -PSEM at ~330 s
DataBlock(28).FfPulseConcat = [RootAnalysisDir,'2016-06-08-0/data003-5/data003-5'] ; % 
DataBlock(28).FfPulseConcatDrugTimes = [360,1730] ; % 
DataBlock(28).DsPath = [RootAnalysisDir,'2016-06-08-0/data001/data001'] ; % moving bars in multidirections for finding DS cells
DataBlock(28).mapEiFlag = true ;

% Control animal, CAUTION: likely Wt
% 50/50 mirror with 0.3 NDF in front (attenuation ~2)
% seemed to stable 
% The Ondansetron condition also slightly changed osmolarity
DataBlock(29).BwPath{1} = [RootAnalysisDir,'2016-06-22-0/data000/data000'] ; % control white noise (BW 15-2) 
DataBlock(29).RgPath{1} = [RootAnalysisDir,'2016-06-22-0/data002/data002'] ; % control 
DataBlock(29).RgPath{2} = [RootAnalysisDir,'2016-06-22-0/data004/data004'] ; % +Ondansetron (5HT3 antagonist) solution 
DataBlock(29).RgPath{3} = [RootAnalysisDir,'2016-06-22-0/data006/data006'] ; % +PSEM (5HT3 antagonist) solution
DataBlock(29).RgPath{4} = [RootAnalysisDir,'2016-06-22-0/data008/data008'] ; % wash 
%DataBlock(29).RgConcat = [RootAnalysisDir,'2016-06-22-0/data002-4-8/data002-4-8'] ; % Ondansetron
%DataBlock(29).RgConcat =[RootAnalysisDir,'2016-06-22-0/data002-6-8/data002-6-8'] ; % PSEM
%DataBlock(29).RgConcat = [RootAnalysisDir,'2016-06-22-0/data002-4-6/data002-4-6'] ; % Control,Ondansetron,PSEM
DataBlock(29).RgConcat = [RootAnalysisDir,'2016-06-22-0/data002-4-6-8/data002-4-6-8'] ; % Control,Ondansetron,PSEM
DataBlock(29).FfPulse{1} = [RootAnalysisDir,'2016-06-22-0/data003/data003'] ; %  
DataBlock(29).FfPulse{2} = [RootAnalysisDir,'2016-06-22-0/data005/data005'] ; %   
DataBlock(29).FfPulse{3} = [RootAnalysisDir,'2016-06-22-0/data007/data007'] ; %
%DataBlock(29).FfPulseConcat
DataBlock(29).FfPulseConcatDrugTimes = [] ; % 
DataBlock(29).DsPath = [RootAnalysisDir,'2016-06-22-0/data001/data001'] ; % moving bars in multidirections for finding DS cells
DataBlock(29).mapEiFlag = true ;

% Control animal, CAUTION: likely Cx57-iCre +/- (het)
% 50/50 mirror with 0.3 NDF in front (attenuation ~2)
% seemed to stable 
DataBlock(30).BwPath{1} = [RootAnalysisDir,'2016-06-28-0/data000/data000'] ; % control white noise (BW 15-2) 
DataBlock(30).RgPath{1} = [RootAnalysisDir,'2016-06-28-0/data002/data002'] ; % control 
DataBlock(30).RgPath{2} = [RootAnalysisDir,'2016-06-28-0/data004/data004'] ; % +PSEM (200 nM) 
DataBlock(30).RgPath{3} = [RootAnalysisDir,'2016-06-28-0/data006/data006'] ; % +PSEM (1 uM)
DataBlock(30).RgPath{4} = [RootAnalysisDir,'2016-06-28-0/data008/data008'] ; % +PSEM (10 uM) 
DataBlock(30).RgPath{5} = [RootAnalysisDir,'2016-06-28-0/data010/data010'] ; % wash 
DataBlock(30).RgConcat = [RootAnalysisDir,'2016-06-28-0/data002-4-6-8-10/data002-4-6-8-10'] ; % Control,Ondansetron,PSEM
DataBlock(30).DsPath = [RootAnalysisDir,'2016-06-28-0/data001/data001'] ; % moving bars in multidirections for finding DS cells
DataBlock(30).mapEiFlag = true ;

% Cx57-Cre + AAV7m8-2A-hEF1a-flex-PSAM - contrast reversing gratings  (+/-PSEM)
% 50/50 mirror with 0.3 NDF in front (attenuation ~2)
% great transfection 
DataBlock(31).BwPath{1} = [RootAnalysisDir,'2016-07-01-1/data000/data000'] ; % control white noise (BW 15-2) 
DataBlock(31).RgPath{1} = [RootAnalysisDir,'2016-07-01-1/data002/data002'] ; % control 
DataBlock(31).RgPath{2} = [RootAnalysisDir,'2016-07-01-1/data004/data004'] ; % +PSEM (200 nM) 
DataBlock(31).RgPath{3} = [RootAnalysisDir,'2016-07-01-1/data006/data006'] ; % +PSEM (1 uM)
DataBlock(31).RgPath{4} = [RootAnalysisDir,'2016-07-01-1/data008/data008'] ; % +PSEM (10 uM) 
DataBlock(31).RgPath{5} = [RootAnalysisDir,'2016-07-01-1/data010/data010'] ; % wash 
DataBlock(31).RgConcat = [RootAnalysisDir,'2016-07-01-1/data002-4-6-8-10/data002-4-6-8-10'] ; % Control,Ondansetron,PSEM
DataBlock(31).FfPulse{1} = [RootAnalysisDir,'2016-07-01-1/data003/data003'] ; % +PSEM at ~320
DataBlock(31).FfPulse{2} = [RootAnalysisDir,'2016-07-01-1/data005/data005'] ; % -PSEM at ~320 s
DataBlock(31).FfPulseConcat = [RootAnalysisDir,'2016-07-01-1/data003-5/data003-5'] ; % 
DataBlock(31).FfPulseConcatDrugTimes = [320,1720] ; % 
DataBlock(31).DsPath = [RootAnalysisDir,'2016-07-01-1/data001/data001'] ; % moving bars in multidirections for finding DS cells
DataBlock(31).mapEiFlag = true ;

% Cx57-Cre + AAV7m8-2A-hEF1a-flex-PSAM - drifting gratings and BW repeats (+/-PSEM)
% PLP (GAD cofactor 150uM) in Ames solution
% 50/50 mirror with 0.3 NDF in front (attenuation ~2)
% weak transfection
DataBlock(32).BwPath{1} = [RootAnalysisDir,'2016-07-20-0/data000/data000'] ; % control white noise (BW 15-2) 
DataBlock(32).DgPath{1} = [RootAnalysisDir,'2016-07-20-0/data002/data002'] ; % control 
DataBlock(32).DgPath{2} = [RootAnalysisDir,'2016-07-20-0/data004/data004'] ; % +PSEM (500 nM) 
DataBlock(32).DgPath{3} = [RootAnalysisDir,'2016-07-20-0/data006/data006'] ; % wash
%DataBlock(32).DgPath{4} = [RootAnalysisDir,'2016-07-20-0/data008/data008'] ; % +PSEM (500 nM) 
DataBlock(32).DgConcat = [RootAnalysisDir,'2016-07-20-0/data002-4-6/data002-4-6'] ; % Control,PSEM, wash, PSEM
DataBlock(32).BwRepeat{1} = [RootAnalysisDir,'2016-07-20-0/data003/data003'] ; % +PSEM at ~410s
DataBlock(32).BwRepeat{2} = [RootAnalysisDir,'2016-07-20-0/data005/data005'] ; % -PSEM at ~ 420s
%DataBlock(32).BwRepeat{3} = [RootAnalysisDir,'2016-07-20-0/data007/data007'] ; % +PSEM at ~ 420s
DataBlock(32).BwRepeatConcat = [RootAnalysisDir,'2016-07-20-0/data003-5/data003-5'] ; % 
DataBlock(32).BwRepeatConcatDrugTimes = [410,1620] ; % 
DataBlock(32).DsPath = [RootAnalysisDir,'2016-07-20-0/data001/data001'] ; % moving bars in multidirections for finding DS cells
DataBlock(32).mapEiFlag = true ;


% Cx57-Cre + AAV7m8-2A-hEF1a-flex-PSAM - drifting gratings and BW repeats (+/-PSEM)
% PLP (GAD cofactor 150uM) in Ames solution
% 50/50 mirror with 0.3 NDF in front (attenuation ~2)
% good transfection
DataBlock(33).BwPath{1} = [RootAnalysisDir,'2016-07-22-0/data000/data000'] ; % control white noise (BW 15-2) 
DataBlock(33).DgPath{1} = [RootAnalysisDir,'2016-07-22-0/data002/data002'] ; % control 
DataBlock(33).DgPath{2} = [RootAnalysisDir,'2016-07-22-0/data004/data004'] ; % +PSEM (500 nM) 
DataBlock(33).DgPath{3} = [RootAnalysisDir,'2016-07-22-0/data006/data006'] ; % wash
DataBlock(33).DgConcat = [RootAnalysisDir,'2016-07-22-0/data002-4-6/data002-4-6'] ; % Control,PSEM, wash, PSEM
DataBlock(33).BwRepeat{1} = [RootAnalysisDir,'2016-07-22-0/data003/data003'] ; % +PSEM at ~420s
DataBlock(33).BwRepeat{2} = [RootAnalysisDir,'2016-07-22-0/data005/data005'] ; % -PSEM at ~ 450s
DataBlock(33).BwRepeatConcat = [RootAnalysisDir,'2016-07-22-0/data003-5/data003-5'] ; % 
DataBlock(33).BwRepeatConcatDrugTimes = [420,1650] ; % 
DataBlock(33).DsPath = [RootAnalysisDir,'2016-07-22-0/data001/data001'] ; % moving bars in multidirections for finding DS cells
DataBlock(33).mapEiFlag = true ;

% Cx57-Cre + AAV7m8-IRES-hEF1a-flex-PSAM - drifting gratings and BW repeats (+/-PSEM)
% PLP (GAD cofactor 150uM) in Ames solution
% 50/50 mirror with 0.3 NDF in front (attenuation ~2)
% sparse transfection
DataBlock(34).BwPath{1} = [RootAnalysisDir,'2016-09-07-0/data000/data000'] ; % control white noise (BW 15-2) 
DataBlock(34).DgPath{1} = [RootAnalysisDir,'2016-09-07-0/data002/data002'] ; % control 
DataBlock(34).DgPath{2} = [RootAnalysisDir,'2016-09-07-0/data004/data004'] ; % +PSEM (500 nM) 
DataBlock(34).DgPath{3} = [RootAnalysisDir,'2016-09-07-0/data006/data006'] ; % wash
DataBlock(34).DgConcat = [RootAnalysisDir,'2016-09-07-0/data002-4-6/data002-4-6'] ; % Control,PSEM, wash, PSEM
DataBlock(34).BwRepeat{1} = [RootAnalysisDir,'2016-09-07-0/data003/data003'] ; % +PSEM at ~420s
DataBlock(34).BwRepeat{2} = [RootAnalysisDir,'2016-09-07-0/data005/data005'] ; % -PSEM at ~ 420s
DataBlock(34).BwRepeatConcat = [RootAnalysisDir,'2016-09-07-0/data003-5/data003-5'] ; % 
DataBlock(34).BwRepeatConcatDrugTimes = [420,1620] ; % 
DataBlock(34).DsPath = [RootAnalysisDir,'2016-09-07-0/data001/data001'] ; % moving bars in multidirections for finding DS cells
DataBlock(34).mapEiFlag = true ;

% Control mouse with DG +/-APV (50uM)
% Caution: seems to be a big change in spikes between BWrepeat blocks - prob not drug related
DataBlock(35).BwPath{1} = [RootAnalysisDir,'2016-09-15-0/data000/data000'] ; % control white noise (BW 15-2) 
DataBlock(35).DgPath{1} = [RootAnalysisDir,'2016-09-15-0/data002/data002'] ; % control 
DataBlock(35).DgPath{2} = [RootAnalysisDir,'2016-09-15-0/data004/data004'] ; % +APV (50 uM) 
DataBlock(35).DgPath{3} = [RootAnalysisDir,'2016-09-15-0/data006/data006'] ; % wash
DataBlock(35).DgConcat = [RootAnalysisDir,'2016-09-15-0/data002-4-6/data002-4-6'] ; % Control,APV, wash,
DataBlock(35).BwRepeat{1} = [RootAnalysisDir,'2016-09-15-0/data003/data003'] ; % +APV at ~430s
DataBlock(35).BwRepeat{2} = [RootAnalysisDir,'2016-09-15-0/data005/data005'] ; % -APV at ~ 400s
DataBlock(35).BwRepeatConcat = [RootAnalysisDir,'2016-09-15-0/data003-5/data003-5'] ; % 
DataBlock(35).BwRepeatConcatDrugTimes = [430,1400] ; % 
DataBlock(35).DsPath = [RootAnalysisDir,'2016-09-15-0/data001/data001'] ; % moving bars in multidirections for finding DS cells
DataBlock(35).mapEiFlag = true ;

% Control mouse (effectively c57) +/- YM90K-DART
DataBlock(36).BwPath{1} = [RootAnalysisDir,'2016-09-28-0/data000/data000'] ; % control white noise (BW 15-2) 
DataBlock(36).DgPath{1} = [RootAnalysisDir,'2016-09-28-0/data003/data003'] ; % control 
DataBlock(36).DgPath{2} = [RootAnalysisDir,'2016-09-28-0/data005/data005'] ; % +YM90K-DART (1uM wash period) 
DataBlock(36).DgPath{3} = [RootAnalysisDir,'2016-09-28-0/data007/data007'] ; % +YM90K-DART (later wash period)
DataBlock(36).DgConcat = [RootAnalysisDir,'2016-09-28-0/data003-5-7/data003-5-7'] ; % Control,DART, DART
DataBlock(36).RgPath{1} = [RootAnalysisDir,'2016-09-28-0/data002/data002'] ; % control 
DataBlock(36).RgPath{2} = [RootAnalysisDir,'2016-09-28-0/data006/data006'] ; % +YM90K-DART (1uM wash period) 
DataBlock(36).RgConcat = [RootAnalysisDir,'2016-09-28-0/data002-6/data002-6'] ; % Control,DART
DataBlock(36).BwRepeat{1} = [RootAnalysisDir,'2016-09-28-0/data004/data004'] ; % 
DataBlock(36).BwRepeat{2} = [RootAnalysisDir,'2016-09-28-0/data004/data004'] ; %
DataBlock(36).BwRepeatConcat = [RootAnalysisDir,'2016-09-28-0/data00/data00'] ; % 
DataBlock(36).BwRepeatConcatDrugTimes = [430,1400] ; % 
DataBlock(36).DsPath = [RootAnalysisDir,'2016-09-28-0/data001/data001'] ; % moving bars in multidirections for finding DS cells
DataBlock(36).mapEiFlag = true ;

% Cx57-Cre + AAV7m8-IRES-hEF1a-flex-PSAM - drifting gratings and BW repeats (+/-PSEM)
% 50/50 mirror with 0.3 NDF in front (attenuation ~2)
% very good transfection
DataBlock(37).BwPath{1} = [RootAnalysisDir,'2016-10-12-0/data000/data000'] ; % control white noise (BW 15-2) 
DataBlock(37).DgPath{1} = [RootAnalysisDir,'2016-10-12-0/data002/data002'] ; % control 
DataBlock(37).DgPath{2} = [RootAnalysisDir,'2016-10-12-0/data004/data004'] ; % +PSEM (500 nM) 
DataBlock(37).DgPath{3} = [RootAnalysisDir,'2016-10-12-0/data006/data006'] ; % wash
DataBlock(37).DgConcat = [RootAnalysisDir,'2016-10-12-0/data002-4-6/data002-4-6'] ; % Control,PSEM, wash, PSEM
DataBlock(37).BwRepeat{1} = [RootAnalysisDir,'2016-10-12-0/data003/data003'] ; % +PSEM at ~420s
DataBlock(37).BwRepeat{2} = [RootAnalysisDir,'2016-10-12-0/data005/data005'] ; % -PSEM at ~ 420s
DataBlock(37).BwRepeatConcat = [RootAnalysisDir,'2016-10-12-0/data003-5/data003-5'] ; % 
DataBlock(37).BwRepeatConcatDrugTimes = [420,1620] ; % 
DataBlock(37).DsPath = [RootAnalysisDir,'2016-10-12-0/data001/data001'] ; % moving bars in multidirections for finding DS cells
DataBlock(37).mapEiFlag = true ;

% c57 control mouse for DART (no DART ever used in this experiment - just simulated protocol)
% using Photons for stimuli!
% 50/50 mirror with 0.3 NDF in front (attenuation ~2)
DataBlock(38).BwPath{1} = [RootAnalysisDir,'2016-11-09-0/data000/data000'] ; % control white noise (BW 15-2) 
DataBlock(38).DgPath{1} = [RootAnalysisDir,'2016-11-09-0/data003/data003'] ; % control 
DataBlock(38).DgPath{2} = [RootAnalysisDir,'2016-11-09-0/data006/data006'] ; % +'DART' (post 'DART', wash) 
DataBlock(38).DgConcat = [RootAnalysisDir,'2016-11-09-0/data003-6/data003-6'] ; % Control,wash
DataBlock(38).RgPath{1} = [RootAnalysisDir,'2016-11-09-0/data002/data002'] ; % control 
DataBlock(38).RgPath{2} = [RootAnalysisDir,'2016-11-09-0/data007/data007'] ; % + 'DART' (post 'DART', wash) 
DataBlock(38).RgConcat = [RootAnalysisDir,'2016-11-09-0/data002-7/data002-7'] ; % Control,DART
DataBlock(38).BwRepeat{1} = [RootAnalysisDir,'2016-11-09-0/data004-5/data004-5'] ; % 
DataBlock(38).BwRepeatConcat = [RootAnalysisDir,'2016-11-09-0/data003-5/data003-5'] ; % 
DataBlock(38).BwRepeatConcatDrugTimes = [420,1280] ; % 
DataBlock(38).DsPath = [RootAnalysisDir,'2016-11-09-0/data001/data001'] ; % moving bars in multidirections for finding DS cells
DataBlock(38).mapEiFlag = true ;

% cx57-Cre HaloTag-injected mouse with DART 
% Good tranfection
% clear rundown later in recording
% 50/50 mirror with 0.3 NDF in front (attenuation ~2)
DataBlock(39).BwPath{1} = [RootAnalysisDir,'2017-07-26-0/data001_JC/data001_JC'] ; % control white noise (BW 15-2) 
DataBlock(39).BwPath{2} = [RootAnalysisDir,'2017-07-26-0/data007/data007'] ; % +DART (BW 15-2) 
DataBlock(39).DgPath{1} = [RootAnalysisDir,'2017-07-26-0/data003/data003'] ; % control 
DataBlock(39).DgPath{2} = [RootAnalysisDir,'2017-07-26-0/data005/data005'] ; % +'DART' (post 'DART', wash) 
DataBlock(39).DgConcat = [RootAnalysisDir,'2017-07-26-0/data003-5/data003-5'] ; % Control,wash
DataBlock(39).RgPath{1} = [RootAnalysisDir,'2017-07-26-0/data002/data002'] ; % control 
DataBlock(39).RgPath{2} = [RootAnalysisDir,'2017-07-26-0/data006/data006'] ; % + 'DART' (post 'DART', wash) 
DataBlock(39).RgConcat = [RootAnalysisDir,'2017-07-26-0/data002-6/data002-6'] ; % Control,DART
DataBlock(39).BwRepeat{1} = [RootAnalysisDir,'2017-07-26-0/data004/data004'] ; % control, +DART, wash
%DataBlock(39).BwRepeat{1} =[RootAnalysisDir,'2017-07-26-0/data009/data009'] ; % + DART
DataBlock(39).BwRepeatConcat = [RootAnalysisDir,'2017-07-26-0/data004/data004'] ; % 
DataBlock(39).BwRepeatConcatDrugTimes = [422,1310] ; % 
DataBlock(39).DsPath{1} = [RootAnalysisDir,'2017-07-26-0/data000/data000'] ; % moving bars in multidirections for finding DS cells
DataBlock(39).DsPath{2} = [RootAnalysisDir,'2017-07-26-0/data008/data008'] ; % + DART
DataBlock(39).mapEiFlag = true ;

% cx57-Cre HaloTag-injected mouse with DART 
% No clear tranfection
% 50/50 mirror with O NDF in front (NO attenuation UNLIKE PREVIOUS EXPERIMENTS)***
DataBlock(40).BwPath{1} = [RootAnalysisDir,'2017-10-06-0/data001/data001'] ; % control white noise (BW 10-2) 
DataBlock(40).BwPath{2} = [RootAnalysisDir,'2017-10-06-0/data007/data007'] ; % +DART (BW 15-2) 
DataBlock(40).DgPath{1} = [RootAnalysisDir,'2017-10-06-0/data003/data003'] ; % control 
DataBlock(40).DgPath{2} = [RootAnalysisDir,'2017-10-06-0/data005/data005'] ; % +'DART' (post 'DART', wash) 
DataBlock(40).DgConcat = [RootAnalysisDir,'2017-10-06-0/data003-5/data003-5'] ; % Control,wash
DataBlock(40).RgPath{1} = [RootAnalysisDir,'2017-10-06-0/data002/data002'] ; % control 
DataBlock(40).RgPath{2} = [RootAnalysisDir,'2017-10-06-0/data006/data006'] ; % + 'DART' (post 'DART', wash) 
DataBlock(40).RgConcat = [RootAnalysisDir,'2017-10-06-0/data002-6/data002-6'] ; % Control,DART
DataBlock(40).BwRepeat{1} = [RootAnalysisDir,'2017-10-06-0/data004/data004'] ; % control, +DART, wash
%DataBlock(40).BwRepeat{1} =[RootAnalysisDir,'2017-10-06-0/data009/data009'] ; % + DART
DataBlock(40).BwRepeatConcat = [RootAnalysisDir,'2017-10-06-0/data004/data004'] ; % 
DataBlock(40).BwRepeatConcatDrugTimes = [422,1310] ; % 
DataBlock(40).DsPath{1} = [RootAnalysisDir,'2017-10-06-0/data001/data001'] ; % moving bars in multidirections for finding DS cells
DataBlock(40).mapEiFlag = true ;

% pv-Cre HaloTag-injected mouse with DART (Genetics show it was Ai14+!!!!)
% Very strong and well spread flourescents but was Ai14+ so unclear if there was transfection
% 50/50 mirror with O NDF in front (NO attenuation UNLIKE PREVIOUS EXPERIMENTS)***
DataBlock(41).BwPath{1} = [RootAnalysisDir,'2017-10-27-0/data001/data001'] ; % control white noise (BW 20-2) 
DataBlock(41).BwPath{2} = [RootAnalysisDir,'2017-10-27-0/data007/data003'] ; % +DART (BW 20-2) 
DataBlock(41).BwRepeat{1} = [RootAnalysisDir,'2017-10-27-0/data002/data002'] ; % control, +DART, wash (+DART)
DataBlock(41).BwRepeat{2} = [RootAnalysisDir,'2017-10-27-0/data005/data005'] ; % control (post +DART), different Ames beaker, back to first
DataBlock(41).BwRepeatConcatDrugTimes = [560,1500] ; % 
%DataBlock(41).BwRepeatConcatDrugTimes = [480,1390] ; % for second BW repeat 
DataBlock(41).DsPath{1} = [RootAnalysisDir,'2017-10-27-0/data000/data000'] ; % control
DataBlock(41).DsPath{2} = [RootAnalysisDir,'2017-10-27-0/data004/data004'] ; % + DART
DataBlock(41).mapEiFlag = true ;

% cx57-Cre HaloTag injected from UCLA
% clear tranfection but little over the array
DataBlock(42).BwPath{1} = [RootAnalysisDir,'2017-11-08-0/data002/data002']  ; % lots more data - check notes

% cx57-Cre HaloTag injected from UCLA
% little transfection over the array but may have been bleached
% lost recording during drug wash in - suspect bubble
DataBlock(43).BwPath{1} = [RootAnalysisDir,'2017-11-09-0/data000/data000']  ; % lots more data - check notes

% cx57-Cre HaloTag injected from UCLA
% little transfection over the array but may have been bleached
DataBlock(44).BwPath{1} = [RootAnalysisDir,'2017-11-15-0/data002/data002']  ; % lots more data - check notes

% Chat-Cre HaloTag injected
% entire array covered in td-tom+ cells 
DataBlock(45).BwPath{1} = [RootAnalysisDir,'2017-11-30-0/data000/data000']  ; % BW 10-2
DataBlock(45).BwPath{2} = [RootAnalysisDir,'2017-11-30-0/data004/data004']  ; % BW 10-2 + DART
DataBlock(45).DsPath{1} = [RootAnalysisDir,'2017-11-30-0/data001/data001'] ; % control
DataBlock(45).DsPath{2} = [RootAnalysisDir,'2017-11-30-0/data003/data003'] ; % + DART
DataBlock(45).DsPathConcat{1} = [RootAnalysisDir,'2017-11-30-0/data001-3/data001-3'] ; % control - +DART
DataBlock(45).BwRepeat{1} = [RootAnalysisDir,'2017-11-30-0/data002/data002'] ; % control (post +DART), different Ames beaker, back to first
DataBlock(45).BwRepeatConcatDrugTimes = [644,1645] ; % 

% PV-Cre HaloTag injected
% little to no transfection over the array
DataBlock(46).BwPath{1} =  [RootAnalysisDir,'2017-12-19-0/data000/data000']  ; % BW 20-2
DataBlock(46).BwPath{2} =  [RootAnalysisDir,'2017-12-19-0/data002/data002']  ; % BW 20-2 (+DART)
DataBlock(46).MbPath{1} =  [RootAnalysisDir,'2017-12-19-0/data001/data001']  ; % moving bar (-/+ DART)
DataBlock(46).MbPathDrugTimes{1} =  [] ; %(+DART)

% Chat-Cre HaloTag injected
% entire array covered in td-tom+ cells
% very sensitive recording
% used Hex to block Ach, DART to block AmpaR in SAC, UBP310 to block Kainate
DataBlock(47).BwPath{1} =  [RootAnalysisDir,'2018-01-10-0/data000/data000'] ; % BW 10-2 (control)
DataBlock(47).DsPath{1} = [RootAnalysisDir,'2018-01-10-0/data001/data001'] ;  % control
DataBlock(47).DsPath{2} = [RootAnalysisDir,'2018-01-10-0/data003/data003'] ;  % control (+Hex, +DART, +UBP310)
DataBlock(47).MbPath{1} =  [RootAnalysisDir,'2018-01-10-0/data002/data002'] ; % moving bar (-/+ Hex, DART, UBP310)
DataBlock(47).MbPathDrugTimes{1} =  [790,2000,2900,4706] ; %(control,+Hex,+DART (still Hex),DART wash (still Hex),+UBP310 (still Hex))
DataBlock(47).PutativeDs = [6,10,17,24,27,30,32,41,47,50,51,57,60,62,68,77,...
    80,88,94,106,112,113,120,123,142,153,158,159,160, 164,186,189,190,194,...
    196,202,209,225,228,235,240] ; % putative DS from qualitative observation of mb

% wt C57 mouse no injection
% Ok recording - spont rate was very high but response was strong
DataBlock(48).BwPath{1} =  [RootAnalysisDir,'2018-01-26-0/data000/data000'] ; % BW 10-2 (control)
DataBlock(48).DsPath{1} = [RootAnalysisDir,'2018-01-26-0/data001/data001'] ;  % control
DataBlock(48).DsPath{2} = [RootAnalysisDir,'2018-01-26-0/data010/data010'] ;  % wash (post UBP 310 wash)
DataBlock(48).MbPath{1} =  [RootAnalysisDir,'2018-01-26-0/data009/data009'] ; % moving bar
DataBlock(48).MbPathDrugTimes{1} =  [950,1916] ; %(control,+UBP310 (Kainate blocker), wash)
DataBlock(48).PutativeDs = [12,21,39,40,41,50,51,101,109,112,163,167,199,209,226,242,244,273] ; % putative DS from qualitative observation of mb

% Chat-Cre HaloTag injected
% ~1/3 array covered in td-tom+ cells
% good recording
% DART to block AmpaR in SAC, then GYKI 52466 to block all Ampa R, then GYKI + UBP310 to block all Ampa and Kainate
DataBlock(49).BwPath{1} =  [RootAnalysisDir,'2018-02-08-0/data000/data000'] ; % BW 10-2 (control)
DataBlock(49).DsPath{1} = [RootAnalysisDir,'2018-02-08-0/data001/data001'] ;  % control
DataBlock(49).MbPath{1} =  [RootAnalysisDir,'2018-02-08-0/data002/data002'] ; % moving bar (+/- DART,GYKI,and UBP310)
%DataBlock(49).MbPath{2} =  [RootAnalysisDir,'2018-02-08-0/data003/data003'] ; % moving bar (start trial shortly after begining drug wash)
DataBlock(49).MbPathDrugTimes{1} =  [830,1750,2900,4500,5010] ; %(control, +DART, wash, +GYKI, wash, +GYKI/UBP310)
DataBlock(49).PutativeDs = [16,87,88,93,94,95,99,100,101,102,103,105,107,109,...
    112,126,163,165,179,201,223,234,263,273,288,306,307,308,312,313,317] ; % putative DS from qualitative observation of mb
DataBlock(49).EpiImagePath = [RootImageDir, '2018-02-08-0/10x tdtom 1.jpg'] ;

% C57/bl6
% ok recording
DataBlock(50).DsPath{1} = [RootAnalysisDir,'2018-05-04-0/data000/data000'] ;  % control
DataBlock(50).MbPath{1} =  [RootAnalysisDir,'2018-05-04-0/data002/data002'] ; % moving bar (+/- GYKI)
DataBlock(50).MbPathDrugTimes{1} =  [1750,3630] ; %(control, +GYKI, wash)
DataBlock(50).PutativeDs = [10,40,49,64,69,84,88,109,110,113,120,124,126,...
    127,133,137,157,164,166,167,184,194,224,227,237,242,268,285,293,...
    297,299,315,316,321,331,348,357,360,361,363,369,377,395] ; % putative DS from qualitative observation of mb

% Chat-Cre HaloTag injected
% entire array covered in td-tom+ cells
% used Hex to block Ach, YM90k-DART(2.0) to block AmpaR in SAC, UBP310 to block Kainate
DataBlock(51).BwPath{1} = [RootAnalysisDir,'2018-10-03-0/data000/data000'] ;  % control
DataBlock(51).MbPath{1} =  [RootAnalysisDir,'2018-10-03-0/data001/data001'] ; % moving bar (control, +Hex, +Hex/DART, +hex/wash, +UBP310)
DataBlock(51).MbPathDrugTimes{1} =  [780,2550,3590,5210] ; %(control, +Hex, +Hex/DART, +hex/wash, +UBP310) +control at 5960
DataBlock(51).PutativeDs = [8,14,15,27,34,35,37,39,50,52,54,57,61,63,65,66,73,86,107,...
    117,120,128,136,161,162,201,204,225,227,233,234,235,243,251,268,269,284,286,292,...
    295,301,308,315,316,358,362,372] ;

% Chat-Cre HaloTag injected
% No tdTom+ cells visible in recorded piece (or anywhere in eyes)
% Ran sham drug (used Ames)
% Mb had multiple speeds, contrasts, and directions
% missing a trigger - causes problems aligning spikes for last several repeats
DataBlock(52).MbPath{1} =  [RootAnalysisDir,'2018-10-25-0/data000/data000'] ; % moving bar (control, +Ames (sham drug), wash)
DataBlock(52).MbPathDrugTimes{1} =  [1910,2790] ; % (control, sham, wash)
%DataBlock(52).PutativeDs =

% Chat-Cre HaloTag injected
% sparse (~20) tdTom+ cells over array
% Ran Gabazine-Dart 2.0
% Mb had multiple speeds, 1 contrast, and directions
DataBlock(53).MbPath{1} =  [RootAnalysisDir,'2018-11-09-0/data000/data000'] ; % moving bar (control, +Dart, wash)
DataBlock(53).MbPathDrugTimes{1} =  [1044,2001] ; % (control, sham, wash)
DataBlock(53).PutativeDs = [14,17,35,37,56,57,63,80,111,113,145,203,216,217,218,221,...
    222,230,233,234,235,241,255,270,275,276,277,281,313,314,324,325,327,341,342,363,391]

% PV-Cre HaloTag injected
% sparse transfection over the array
% poor recording
% Ran Gabazine-Dart 2.0
DataBlock(54).BwPath{1} = [RootAnalysisDir,'2018-11-14-0/data000/data000'] ; % pre Dart
DataBlock(54).BwPath{2} = [RootAnalysisDir,'2018-11-14-0/data004/data004'] ; % +Dart
DataBlock(54).DgPath{1} = [RootAnalysisDir,'2018-11-14-0/data001/data001'] ; % pre Dart
DataBlock(54).DgPath{2} = [RootAnalysisDir,'2018-11-14-0/data003/data003'] ; % +Dart
DataBlock(54).BwRepeat{1} = [RootAnalysisDir,'2018-11-14-0/data002/data002'] ; % control, Dart, Dart wash 
DataBlock(54).BwRepeatDrugTimes = [970,1845] ; %
DataBlock(54).EpiImagePath = [RootImageDir,'2018-11-14-0/10x td-tom 6.jpg'] ;

% PCP2-Cre HT injected
% YM90K-DART
DataBlock(55).BwPath{1} = [RootAnalysisDir,'2018-11-16-0/data002/data002'] ; % control,
DataBlock(55).BwPath{2} = [RootAnalysisDir,'2018-11-16-0/data004/data004'] ; % post DART
DataBlock(55).DgPath{1} = [RootAnalysisDir,'2018-11-16-0/data001/data001'] ; % pre Dart
DataBlock(55).DgPath{2} = [RootAnalysisDir,'2018-11-16-0/data005/data005'] ; % +Dart
DataBlock(55).FfPulse{1} = [RootAnalysisDir,'2018-11-16-0/data001/data000'] ; % pre Dart
DataBlock(55).FfPulse{2} = [RootAnalysisDir,'2018-11-16-0/data005/data006'] ; % +Dart
DataBlock(55).BwRepeat{1} = [RootAnalysisDir,'2018-11-16-0/data003/data003'] ; % control, Dart, Dart wash 
DataBlock(55).BwRepeatDrugTimes = [808,1580] ; %
DataBlock(55).EpiImagePath = [RootImageDir,'2018-11-16-0/td-tom 4x.jpg'] ;

% PV-Cre HT injected
% Gabazine+Cy5 - DART
% ok transfection
DataBlock(56).BwPath{1} = [RootAnalysisDir,'2018-11-28-0/data000/data000'] ; % control,BW 20-2
DataBlock(56).BwPath{2} = [RootAnalysisDir,'2018-11-28-0/data004/data004'] ; % post DART,BW 20-2
DataBlock(56).DgPath{1} = [RootAnalysisDir,'2018-11-28-0/data001/data001'] ; % pre Dart
DataBlock(56).DgPath{2} = [RootAnalysisDir,'2018-11-28-0/data003/data003'] ; % +Dart
DataBlock(56).DgConcat{1} = [RootAnalysisDir,'2018-11-28-0/data001-3/data001-3'] ; % control/+Dart
DataBlock(56).BwRepeat{1} = [RootAnalysisDir,'2018-11-28-0/data002/data002'] ; % control, Dart, Dart wash 
DataBlock(56).BwRepeatDrugTimes = [958,1910] ; %
DataBlock(56).IrImagePath = [RootImageDir,'2018-11-28-0/10x tdTom 1.jpg'] ;
DataBlock(56).EpiImagePath = [RootImageDir,'2018-11-28-0/10x tdTom 1.jpg'] ;

% PV-Cre HT injected
% Gabazine+Cy5 - DART
% ok transfection
DataBlock(57).BwPath{1} = [RootAnalysisDir,'2018-12-05-0/data000/data000'] ; % control,BW 10-2
DataBlock(57).BwPath{2} = [RootAnalysisDir,'2018-12-05-0/data006/data006'] ; % post DART,BW 10-2
DataBlock(57).DgPath{1} = [RootAnalysisDir,'2018-12-05-0/data002/data002'] ; % pre Dart
DataBlock(57).DgPath{2} = [RootAnalysisDir,'2018-12-05-0/data004/data004'] ; % +Dart
DataBlock(57).DgConcat{1} = [RootAnalysisDir,'2018-12-05-0/data002-4/data002-4'] ; % control/+Dart
DataBlock(57).BwRepeat{1} = [RootAnalysisDir,'2018-12-05-0/data003/data003'] ; % control, Dart, Dart wash 
DataBlock(57).BwRepeatDrugTimes = [960,1912] ; % (Dart,wash)
DataBlock(57).FfPulse{1} = [RootAnalysisDir,'2018-12-05-0/data001/data000'] ; % pre Dart
DataBlock(57).FfPulse{2} = [RootAnalysisDir,'2018-12-05-0/data005/data006'] ; % +Dart
DataBlock(57).IrImagePath = [RootImageDir,'2018-12-05-0/10x RFP 5.jpg'] ;
DataBlock(57).EpiImagePath = [RootImageDir,'2018-12-05-0/10x RFP 5.jpg'] ;

% Params defaults
Params.empty = [] ;
ForIgor = [] ;

% run functions

if perform.KoBwAnalyzer2 ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = KoBwAnalyzer2(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.KoBwAnalyzer3 ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = KoBwAnalyzer3(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.KoBwAnalyzer4 ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = KoBwAnalyzer4(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end


if perform.KoStepAnalyzer ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = KoStepAnalyzer(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.KoDriftingGratingAnalyzer ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = KoDriftingGratingAnalyzer(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end
 
if perform.KoSpatialTuningAnalyzer ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = KoSpatialTuningAnalyzer(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.KoSpatialTuningAnalyzerNoMapping ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = KoSpatialTuningAnalyzerNoMapping(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.KoSpatialTuningAnalyzerNoMappingV2 ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = KoSpatialTuningAnalyzerNoMappingV2(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.KoSpatialTuningAnalyzerConcat ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = KoSpatialTuningAnalyzerConcat(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.KoSpatialTuningAnalyzerConcatPlusMapping ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = KoSpatialTuningAnalyzerConcatPlusMapping(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.KoSpatialTuningAnalyzerConcatPlusMappingV2 ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = KoSpatialTuningAnalyzerConcatPlusMappingV2(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.KoSpatialTuningAnalyzerConcatPlusMappingV3 ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = KoSpatialTuningAnalyzerConcatPlusMappingV3(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.KoSpatialTuningAnalyzerConcatPlusMappingV4 ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) 
        Temp = KoSpatialTuningAnalyzerConcatPlusMappingV4(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end


if perform.KoStepConcatAnalyzer ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = KoStepConcatAnalyzer(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.KoStepConcatAnalyzerPlusMapping ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = KoStepConcatAnalyzerPlusMapping(DataBlock, DB, Params) ;
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


if perform.KoReverseGratingsAnalyzerConcatPlusMappingV1 ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = KoReverseGratingsAnalyzerConcatPlusMappingV1(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.KoReverseGratingsAnalyzerConcatPlusMappingMultiSetV1 ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = KoReverseGratingsAnalyzerConcatPlusMappingMultiSetV1(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.KoGratingsAnalyzerConcatPlusMappingMultiSetV1 ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = KoGratingsAnalyzerConcatPlusMappingMultiSetV1(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.KoBwRepeatConcatAnalyzerPlusMapping ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = KoBwRepeatConcatAnalyzerPlusMapping(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end
    

if perform.rfProfileAnalysis ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = rfProfileAnalysis(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.KoDsConcatAnalysis ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = KoDsConcatAnalysis(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

if perform.MbDrugAnalysis ;
    for a=1:length(DBset{DBset_id}) ; % for each data block in set
        DB = DBset{DBset_id}(a) ;
        Temp = MbDrugAnalysis(DataBlock, DB, Params) ;
        ForIgor = mergeStruct(ForIgor,Temp) ;
    end
end

end % function DataBlocks_KO
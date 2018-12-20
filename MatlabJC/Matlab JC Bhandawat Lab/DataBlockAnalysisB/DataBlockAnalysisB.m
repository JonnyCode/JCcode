%Data notes, data blocks and analysis logicals
% copied from "DataBlockAnalysis" but have changed so data block loop (A) is performed in called functions,
% this necisitates functions to be modified.
%Bhandawat lab
%J Cafaro 1/17/13

%% Analysis logicals (true or false) and data blocks to be analyzed

perform.PNadaptationAnalysisI = false ; % similar to H but using performs data block loop and population analysis within

perform.PNadaptationAnalysisJ = false ; % similar to I but focuses on psth peak, calculates bg and rest differently, and calc discriminability

perform.PNadaptationAnalysisK = false ; % similar to J but added divisive normalization

perform.PNadaptationAnalysisL = false ; % similar to K but changed normalization to subtractive

perform.adaptationAnalysisApn = false ; % similar to K but combines ORN and PN codes

perform.adaptationAnalysisBpn = false ; % similar to A but focuses on normalization factors 

perform.adaptationAnalysisDpn = false ; % refined version of B

perform.adaptationAnalysisEpn = false ; % refined version of D

perform.adaptationAnalysisFpn = true ; % refined version of E


perform.ORNadaptationAnalysisA = false ; % adapted from PNadaptationAnalysisI but for ORNs

perform.ORNadaptationAnalysisB = false ; % adapted from PNadaptationAnalysisJ but for ORNs

perform.ORNadaptationAnalysisC = false ; % adapted from PNadaptationAnalysisK but for ORNs

perform.adaptationAnalysisAorn = false ; % same as adaptationAnalysisA but runs with ORN data flag for spike detection

perform.adaptationAnalysisBorn = false ; % same as adaptationAnalysisB but runs with ORN data flag for spike detection

perform.adaptationAnalysisDorn = false ; % refined version of B

perform.adaptationAnalysisEorn = false ; % refined version of D

perform.adaptationAnalysisForn = false ; % refined version of E


perform.PIDAnalysisA = false ; % analysis of PID recordings for +/- background DB

% perform analysis of following data blocks...
%A = [66,69,72:91,95:112] ; % ORN DB
%A = [46:51,53:58,60:61,92:94,124:131,133:139] ; % PN DB
%A = [113:119,121:123] ; % ORN DB (ethyl acetate)
%A = [140:149] ; % PN DB (ethyl acetate)
%A = [150,151] ; % PID DB (2-butanone)
A=[46] ;
%% Data blocks and notes

%putative PN, 
%in around -45mV, 
%current data will not be correct,
%cage was not down so a lot of 60 hz noise
%other potentially usefull trials
Input(1).cellname = '111229_2' ; 
Input(1).Spont = '[5:34,37:39]' ; % changing steady injected current
 
%putative PN,
%in around -33mV
%current data will not be correct
Input(2).cellname = '111229_3' ; 
Input(2).Spont = '[0:15]' ; % changing steady injected current

%putative PN,
%in around -50mV
%current data will not be correct
%other data playing with strong current injections (trials 36:45)
%other data playing with injecting current ramp 
%other data playing with lights in room and puffing air on fly
Input(3).cellname = '111230_1' ; 
Input(3).Spont = '[5:25]' ; % changing steady injected current

%putative PN,
%in around -35mV, okay seal
%other data playing with floor vibrations
%other data playing with clapping - cell hyperpol during clapping
%ramped current until cell crashed (trial 30)
Input(4).cellname = '111230_2' ; 
Input(4).Spont = '[1:10]' ; 

%putative PN, gfp+
%in around -54mV
%other data playing with clapping
%ramped current until cell crashed (trials 18)
%this cell had big swings in potential
Input(5).cellname = '111230_3' ; 
Input(5).Spont = '[0:10]' ;

% ORN ab4 from Seugh-Hye in OR83b+/-
% assing spontaneous rate and odor response
Input(6).cellname = '120313#1' ;
Input(6).Spont = '[6:9]' ; 
Input(6).OdorRsp = '[10]' ; % trans-2-hexanol at 10^-4

% ORN ab4 from Seugh-Hye in OR83b+/-
% assing spontaneous rate and odor response
Input(7).cellname = '120313#2' ;
Input(7).Spont = '[0:4]' ; 
Input(7).OdorRsp = '[5]' ; % trans-2-hexanol at 10^-4

% ORN ab4 from Seugh-Hye in wt
% assing spontaneous rate and odor response
Input(8).cellname = '120314#3' ;
Input(8).Spont = '[0:3]' ; 
Input(8).OdorRsp = '[5]' ; % trans-2-hexanol at 10^-4

% putative PN with I injection
Input(9).cellname = '040412#2' ;
Input(9).Iinj = '[1:x]' ;

% putative PN with I injection
Input(10).cellname = '040412#3' ;
Input(10).Iinj = '[1:xx]' ;

% VM7 PN with 10^-4 Butatnone pulse +/- background 
% medial, gfp+, strong butanone response and after hyperpolarization
% one antenna was ripped off
Input(11).cellname = '040912#1' ;
Input(11).OdorRsp = '[]' ; % after around epoch 2 of 5 sec trials put on background
Input(11).OdorRspPlusBg = '[]' ; 

% putative PN with 10^-4 butanone pulse and long pulses of 10^-2 ethly
% acetate, medial 
% one antenna was ripped off
Input(12).cellname = '040912#2' ;
Input(12).OdorRsp = '[]' ; % trial 22 ethyl acetate?, epoch 28 no odor?

% VM7 PN with 10^-4 2-butanone pulse +/- background
% gfp+ medial PN with strong response
Input(13).cellname = '120423_1' ; 
Input(13).OdorRsp = '[0:95,97:134]' ; % odor response +/- background
Input(13).OdorBgTimes{1} = {'[23-Apr-2012 14:38:33]','[23-Apr-2012 14:44:48]','[23-Apr-2012 14:52:19]','[23-Apr-2012 14:59:00]'} ; % time at which backgrounds 10^-7 2 butanone were present, and then removed
Input(13).OdorBgTimes{2} = {'[23-Apr-2012 15:13:33]','[23-Apr-2012 15:26:55]'} ; % time at which backgrounds 10^-6 2 butanone were present, and removed

% VM7 PN (gfp+ strong 2b response),
% with 10^-4 2-butanone pulse +/- pure low flow (12 ml/min) background
Input(14).cellname = '120430_1' ;
Input(14).OdorRsp = '[1:13]' ; % odor response +/- background
Input(14).OdorBgTimes{1} = {'[30-Apr-2012 10:31:45]','[30-Apr-2012 10:40:03]'} ; % 10 ml/min
Input(14).OdorBgTimes{2} = {'[30-Apr-2012 10:40:03]'} ; % 2-3 ml min

% puatative PN (not VM7)
% odor valve noise
Input(15).cellname = '120511_1' ;
Input(15).ValveNoise = '[21:25]' ;

% VM7 PN (gfp+ strong 2b response),
% with 10^-4 2-butatnone pulse +/- pure low flow (0, 2.5ml/min) background
Input(16).cellname = '120516_1' ;
Input(16).OdorRsp = '[2:62]' ;
Input(16).OdorBgTimes{1} = {'[16-May-2012 14:03:11]','[16-May-2012 14:31:26]'} ; % no odor
Input(16).OdorBgTimes{2} = {'[16-May-2012 13:42:48]','[16-May-2012 13:56:38]','[16-May-2012 14:10:10]','[16-May-2012 14:24:23]'} ; % pure at around 0 ml/min
Input(16).OdorBgTimes{3} = {'[16-May-2012 13:50:23]','[16-May-2012 14:17:22]'} ; % pure at around 2-3 ml/mi
Input(16).OdorBgTimes{4} = {'[16-May-2012 14:38:21]'} ; % no odor but depolarizing cell

% VM7 PN (gfp+ strong 2b response)
% with 10^-4 2-butanone or 10^-4 ethyl acetate pulse +/- pure low flow 2-butanone (0 ml/min) background
Input(17).cellname = '120518_1' ;
Input(17).OdorRsp = '[8:24]' ; % ethyl acetate pulse, butanone background
%Input(17).OdorRsp = '[25:30]' ; % 2-butanone pulse but few spikes
Input(17).OdorBgTimes{1} = {'[18-May-2012 12:19:30]'} ; % no odor
Input(17).OdorBgTimes{2} = {'[18-May-2012 12:12:43]'} ; % pure at 0ml/min

% VM7 PN (gfp+ strong 2b response)
% with 10^-4 2-butanone pulse +/- pure low flow (0,2ml/min) background
Input(18).cellname = '120518_2' ;
Input(18).OdorRsp = '[0:48]' ;
%Input(18).OdorRsp = '[49:53]' ; % odor pulse - background + hyperpol to see if plataue could be recovered
Input(18).OdorBgTimes{1} = {'[18-May-2012 14:38:25]'} ; % no odor
Input(18).OdorBgTimes{2} = {'[18-May-2012 14:17:35]','[18-May-2012 14:31:35]','[18-May-2012 14:46:48]'} ; % pure at around 0 ml/min
Input(18).OdorBgTimes{3} = {'[18-May-2012 14:24:36]','[18-May-2012 14:55:09]'} ; % pure at 2 ml/mi
Input(18).OdorBgTimes{4} = {'[18-May-2012 15:02:05]'} ; % pure at 3 ml/mi

% ab2 from seugh-hey age 3 Orco-gal4:uas-orco (air 40/200 ml/min, odor on 13-13.5sec)
Input(19).cellname = '120514#1' ;
Input(19).OdorRsp = '[2:6]' ; % ethyl acetate 10^-4
Input(19).OdorRsp = '[2:5]' ; % ethyl acetate 10^-4

%ab2 from seugh-hey age 3 uas-orco (air 20/100 ml/min, odor on 13-13.5sec)
Input(20).cellname = '120514#2' ; 
Input(20).OdorRsp = '[4:8]' ; % ethyl buterate 10^-2

% ab2 from seugh-hey age 3 uas-orco (air 20/100 ml/min, "sensilium didn't look healthy")
Input(21).cellname = '120518#1' ; 
Input(21).OdorRsp = '[9]' ; % ethyl acetate 10^-4 odor pulse 10-11 sec
Input(21).OdorRsp = '[10]' ; % ethyl butyrate 10^-2 odor 13-13.5 sec
Input(21).OdorRsp = '[11:12]' ; % ethyl butyrate 10^-2, 13-14 sec

% ab2? from seugh-hey age 3 orco-ga4; uas-orco (air 20/100 ml/min)
Input(22).cellname = '120518#1' ;
Input(22).OdorRsp = '[20]' ; % ethyl acetate 10^-4, 13-14 sec
Input(22).OdorRsp = '[21:22]' ; % ethyl acetate 10^-4, 13-13.5 sec
%Input(22).OdorRsp = '[23]' ; % ethyl butyrate 10^-2, 13-13.5 sec

% ab2 from seugh-hey age 3 orco-ga4; uas-orco (air 20/100 ml/min)
Input(23).cellname = '120518#1' ;
Input(23).OdorRsp = '[35:39]' ; % ethyl acetate 10^-4, 13-13.5 sec
%Input(23).OdorRsp = '[40:42]' ; % ethyl butyrate 10^-2, 13-13.5 sec

% 120518#1 also has some data on 2x ab3?, ab2?

% puatative PN (not VM7 but gfp+)
% odor valve noise
Input(24).cellname = '120524_1' ;
Input(24).ValveNoise = '[10:11]' ; % ethyl butyrate 10^-2 (different seeds)

% VM7 PN (gfp+ strong 2b subthresh response)
% cell was hyperpol and weak spiking and another cell in the prep was similar 
Input(25).cellname = '120604_1' ;
Input(25).OdorRsp = '[4:36]' ; % 2-butanone multiple concentrations
Input(25).OdorBgTimes{1} = {'[04-Jun-2012 10:58:41]','[04-Jun-2012 11:06:57]','[04-Jun-2012 11:14:29]','[04-Jun-2012 11:21:42]'} ; % no odor but multiple pulse concentrations
Input(25).OdorBgTimes{2} = {'[04-Jun-2012 11:29:02]'} ; % pure 2-butatnone attached 0ml/min
Input(25).OdorBgTimes{3} = {'[04-Jun-2012 11:35:52]'} ; % pure 2-butatnone attached 1 ml/min

% vm7 (gfp? strong 2b response)
% this cell also has some long responses to attachment of pure 2-butanone (99-100)
Input(26).cellname = '120604_1' ;
Input(26).OdorRsp = '[39:98]' ;  % 2-butanone multiple concentrations
Input(26).OdorBgTimes{1} = {'[04-Jun-2012 14:05:18]','[04-Jun-2012 14:12:48]','[04-Jun-2012 14:20:18]','[04-Jun-2012 14:27:36]'} ; % no background
Input(26).OdorBgTimes{2} = {'[04-Jun-2012 14:35:18]','[04-Jun-2012 14:50:48]','[04-Jun-2012 15:00:02]'} ; % pure 2-butatnone attached 0ml/min
Input(26).OdorBgTimes{3} = {'[04-Jun-2012 14:42:15]','[04-Jun-2012 14:49:17]'} ; % pure 2-butatnone attached 1 ml/min
Input(26).OdorBgTimes{4} = {'[04-Jun-2012 15:06:43]','[04-Jun-2012 15:13:33]','[04-Jun-2012 15:14:28]'} ; % pure 2-but <<1ml/min

% vm7 (gpf+ strong 2b response)
% this cell has no background but several concentrations
% timing of the pulses are not always consistent (computer issues)
Input(27).cellname = '120702_1' ;
Input(27).OdorRsp = '[]' ;

% vm7 (gfp+ medium 2b response)
% +/- background at several concentrations
Input(28).cellname = '120716_1' ;
Input(28).OdorRsp{1} = {'[11:16]','[5:10]','[17:22,23:28]','[]'} ; % no background
Input(28).OdorRsp{2} = {'[47:52]','[29:34]','[35:40]','[41:46]'} ; % 5 ml/min large vial at 10^-4 2-but
Input(28).OdorRsp{3} = {'[53:58]','[]','[59:64]','[65:70]'} ; % no background (after background)
Input(28).OdorConcentration = '[10^-5,10^-4,10^-3,10^-2]' ; % 
Input(28).OdorBgTimes = {'[16-Jul-2012 13:34:54]','[16-Jul-2012 14:02:48]'} ;

% vm7 (gfp+ strong 2b response)
% no background at several concentrations
Input(29).cellname = '120802_1' ;
Input(29).OdorRsp{1} = {'[14:19]','[8:13]','[2:7]','[20:25]'} ; % no background
Input(29).OdorConcentration = '[10^-8,10^-7,10^-4,10^-2]' ;
Input(29).OdorBgTimes = {} ;

% Vm7 (gfp+ strong 2b response)
% +/- background at several concentrations (thresh, middle, max)
% some data at higher air rate not listed below
Input(30).cellname = '120813_1' ; 
Input(30).OdorRsp{1} = {'[12:17]','[18:23]','[]','[30:35]','[24:29]'} ; % no background
%Input(30).OdorRsp{1} = {'[12:17]','[18:23]','[]','[6:11,30:35]','[24:29]'} ; % no background
Input(30).OdorRsp{2} = {'[]','[42:47]','[54:59]','[36:41,60:65]','[48:53]'} ; % 10^-4 in big vial at 5ml/min
Input(30).OdorConcentration = '[10^-7,10^-6,10^-5,10^-4,10^-2]' ;
Input(30).OdorBgTimes = {'[13-Aug-2012 11:56:34]'} ; 

% VM7 (weak gfp+, good 2b response)
% +/- background at several concentrations 
Input(31).cellname = '120918_1' ;
Input(31).OdorRsp{1} = {'[2:7,8:13]','[26:31]','[50:55]'} ; % no background
Input(31).OdorRsp{2} = {'[14:19]','[32:37]','[56:61]'} ; % 10^-4 in big vial at 5.5ml/min  
Input(31).OdorRsp{3} = {'[20:25]','[38:43,44:49]','[62:67]'} ; % wash no background
Input(31).OdorConcentration = '[10^-6,10^-5,10^-4]' ;
Input(31).OdorBgTimes = {'[18-Sep-2012 14:15:47]','[18-Sep-2012 14:22:22]','[18-Sep-2012 14:36:06]','[18-Sep-2012 14:42:58]','[18-Sep-2012 15:05:18]','[18-Sep-2012 15:12:20]'} ; 

% VM7 (gfp+ good 2-but response)
% +/- background
% spikes prob not resolvable and cell depolarized during recording
Input(32).cellname = '120924_1' ;
Input(32).OdorRsp{1} = {'[18:23,30:35]'} ; % no background
Input(32).OdorRsp{2} = {'[24:29]'} ; % 10^-4 in big vial at 5.5ml/min    
Input(32).OdorConcentration = '[10^-4]' ;
Input(32).OdorBgTimes = {'[24-Sep-2012 10:47:44]'} ; 

% VM7 (gfp+ good 2-but response)
% +/- background at several concentrations
% very sensetive cell (good 10^-7 response)
Input(33).cellname = '120927_1' ;
Input(33).OdorRsp{1} = {'[5:15]','[48:58]','[81:91]'} ; % no background
Input(33).OdorRsp{2} = {'[16:26]','[59:69]','[92:102]'} ; % 10^-4 in big vial at 5.5ml/min 
Input(33).OdorRsp{3} = {'[27:37]','[70:80]','[103:113]'} ; % no background (WASH)
Input(33).OdorConcentration = '[10^-7,10^-6,10^-2]' ;
Input(33).OdorBgTimes = {'[27-Sep-2012 11:06:32]','[27-Sep-2012 11:18:23]','[27-Sep-2012 11:57:24]','[27-Sep-2012 12:09:49]','[27-Sep-2012 12:40:02]','[27-Sep-2012 12:57:11]'} ; 

% VM7 (gfp+ good 2-but respone)
% +/- background at 1 odor concentration
% cell was depolarized but stable with current injection
% cell was not very sensetive 
Input(34).cellname = '121001_1' ;
Input(34).OdorRsp{1} = {'[12:22,34:44]'} ; % no background
Input(34).OdorRsp{2} = {'[23:33]'} ; % 10^-4 in big vial at 5.5ml/min    
Input(34).OdorConcentration = '[10^-5]' ;
Input(34).OdorBgTimes = {'[01-Oct-2012 11:27:04]','[01-Oct-2012 11:39:41]'} ; 

% VM7 (gfp+ good 2-but response)
% +/- background at several background concentrations
% cell was sensetive (this may be because using fresh odor vials)
Input(35).cellname = '121005_1' ;
Input(35).OdorRsp{1} = {'[3:13]','[36:46]','[69:79]'} ; % no background
Input(35).OdorRsp{2} = {'[14:24]','[47:57]','[80:90]'} ; % 10^-4 in big vial at 5.5ml/min  
Input(35).OdorRsp{3} = {'[25:35]','[58:68]','[91:92]'} ; % no background (WASH)
Input(35).OdorConcentration = '[10^-7,10^-6,10^-4]' ;
Input(35).OdorBgTimes = {'[05-Oct-2012 14:00:57]','[05-Oct-2012 14:13:32]','[05-Oct-2012 14:37:48]','[05-Oct-2012 14:49:23]','[05-Oct-2012 15:14:19]','[05-Oct-2012 15:25:53]'} ; 

% DL5 (gfp+ strong response to trans-2-hex at 10^-2)
Input(36).cellname = '120912_1' ;
Input(36).PulseRsp = '[24:40]' ; % mix of pulse durations and interpulse intervals
Input(36).ValveNoise = '[42:52]' ; % different seeds

% DL5 (gfp+ strong response to trans-2-hex at 10^-7)
Input(37).cellname = '121018_1' ;

% VM7 (gfp+ strong 2-but response)
Input(38).cellname = '121019_1' ;
%Input(38).OdorRsp{1} = {'[4:8]','[41:45]','[61:65]'} ; % no background (put on a background 9-18 but was concerned the odor was old so remade background vial and started again) 
Input(38).OdorRsp{1} = {'[19:25]','[41:45]','[61:65]'} ; % no background
Input(38).OdorRsp{2} = {'[26:35]','[46:55]','[66:75]'} ; % 10^-4 in big vial at 5.5ml/min  
Input(38).OdorRsp{3} = {'[36:40]','[56:60]','[76:77]'} ; % no background (WASH)
Input(38).OdorConcentration = '[10^-7,10^-4,10^-2]' ;
Input(38).OdorBgTimes = {'[19-Oct-2012 10:41:07]','[19-Oct-2012 10:51:55]','[19-Oct-2012 11:00:25]','[19-Oct-2012 11:22:13]','[19-Oct-2012 11:33:08]','[19-Oct-2012 11:44:18]','[19-Oct-2012 11:55:14 ]'} ; 

% VM7 (gfp+ strong 2-but response)
% tested background odors with valve in this cell
Input(39).cellname = '121030_1' ;
Input(39).OdorRsp{1} = '[14:21]' ; % 10^-4 at 5.5ml/min
Input(39).OdorRsp{2} = '[22:25]' ; % 10^-5 at 5.5 ml/min

% VM7 (gfp+ strong 2-but responses)
% multiple concentrations +/- a background 10^-5 at 5.5-6ml/min
Input(40).cellname = '121031_1' ;
Input(40).OdorConcentration = '[10^-7,10^-6,10^-4,10^-3,10^-2]' ;
Input(40).OdorBgTimes = {'[31-Oct-2012 13:29:56]'} ;
Input(40).OdorRsp{1} = {'[4:9]','[10:15]','[16:21]','[22:27]','[28:33]'} ; 
Input(40).OdorRsp{2} = {'[58:63]','[46:51]','[52:57]','[64:69]','[40:45]'} ;

% VM7 (gfp+ strong 2-but response)
% multiple concentrations (no background)
Input(41).cellname = '121101_1' ;

% VM7 (gfp+ strong 2-but response)
% FIRST CELL WITH BACKGROUND ON VALVE
% this cell has high variability in background response caused by exploring conditions impacting odor stimulus
Input(42).cellname = '121114_1' ;

% VM7
% +/- bacgkround on valve at several pulse concentrations
% this cell has some trials at the end exploring residual activity (no carbon filter used)
Input(43).cellname = '121115_1' ;
Input(43).OdorConcentration = '[10^-7,10^-6,10^-4,10^-3,10^-2]' ;
Input(43).OdorRsp{1} = {'[19:24]','[2:7]','[37:42]','[55:60]','[73:78]'} ; % no background
Input(43).OdorRsp{2} = {'[25:30]','[8:13]','[43:48]','[61:66]','[79:84]'} ; % 10^-5 bottle (primed) at 10ml/min
Input(43).OdorRsp{3} = {'[31:36]','[14:18]','[49:54]','[67:72]','[]'} ; % wash

% VM7
% +/- bacgkround on valve at several pulse concentrations
% this cell also has some data exploring carbon filters
Input(44).cellname = '121119_1' ;
Input(44).OdorConcentration = '[10^-7,10^-6,10^-5,10^-4,10^-3]' ;
Input(44).OdorRsp{1} = {'[27:32]','[7:13]','[72:77]','[45:50]','[91:96]'} ; % no background
Input(44).OdorRsp{2} = {'[33:38]','[14:19]','[79:84]','[51:56]','[97:102]'} ; % 10^-5 bottle (primed) at 10ml/min
Input(44).OdorRsp{3} = {'[39:44]','[20:25]','[85:90]','[57:62]','[103:108]'} ; % wash

%VM7
% +/- background on vavle at several pulse concentrations
% this cell has high resting variability
% this cell also has some data with varaible steady current injection
Input(45).cellname = '121121_1' ;
Input(45).OdorConcentration = '[10^-7,10^-6,10^-4,10^-3]' ;
Input(45).OdorRsp{1} = {'[24:31]','[4:10]','[32:37]','[56:61]'} ; % no background
Input(45).OdorRsp{2} = {'[]','[12:17]','[38:43]','[62:67]'} ; % 10^-4 bottle (primed) at 10ml/min
Input(45).OdorRsp{3} = {'[]','[18:23]','[44:55]','[68:73]'} ; % wash
Input(45).OdorRsp{4} = {'[]','[]','[]','[75:80]'} ; % 10^-5 bottle (primed) at 10ml/min

% VM7
% this cell has high resting variability
%*IMPORTANT CHANGES BEYOND THIS DATE*
% +/- background on valve at same flow rate as pulse (100ml/min with constant air at 1L/min) 
% previous work had pulse at 200ml/min and background at variable rates usually around 5-10 ml/min with constant air at 2L/min
% this experiment protocol ('BackgroundOdorAlternate') also alternates background and no background
% charcoal filters were used in these experiments on both background and pulse vials
% this cell used for example in Fig. 1 (6/24/14)
Input(46).cellname = '121203_1' ;
Input(46).OdorConcentration = '[10^-7,10^-6,10^-5,10^-4,10^-3]' ;
Input(46).OdorRsp{1} = {'[11:2:24]','[25:2:42]','[60:2:75]','[43:2:56]','[76:2:91]'} ; % no background
Input(46).OdorRsp{2} = {'[12:2:24]','[26:2:42]','[61:2:75]','[44:2:56]','[77:2:91]'} ; % 10^-7 at 100ml/min
Input(46).BgConcentration = {'0','Log7'} ;

% VM7
% +/- background at several pulse concentrations
Input(47).cellname = '121205_1' 
Input(47).OdorConcentration = '[10^-7,10^-6,10^-5]' ;
Input(47).OdorRsp{1} = {'[4:2:17]','[20:2:33]','[36:2:49]'} ; % no background
Input(47).OdorRsp{2} = {'[5:2:17]','[21:2:33]','[37:2:49]'} ; % 10^-7 at 100ml/min
% Input(47).OdorRsp{1} = {'[4:2:17]','[20:2:33]','[36:2:49]','[52:2:65]','[68:2:73]'} ; % no background
% Input(47).OdorRsp{2} = {'[5:2:17]','[21:2:33]','[37:2:49]','[53:2:65]','[69:2:73]'} ; % 10^-7 at 100ml/min
Input(47).BgConcentration = {'0','Log7'} ;

% vm7
% +/- background at several pulse concentrations
% small spikes
Input(48).cellname = '121206_1' 
Input(48).OdorConcentration = '[10^-6,10^-5,10^-4]' ;
Input(48).OdorRsp{1} = {'[21:2:34]','[35:2:48]','[63:2:76]'} ; % no background
Input(48).OdorRsp{2} = {'[22:2:34]','[36:2:48]','[64:2:76]'} ; % 10^-6 at 100ml/min
Input(48).BgConcentration = {'0','Log6'} ;

% vm7
% +/- background at several pulse concentrations
% small spikes and high rest potential
Input(49).cellname = '121212_1' ;
Input(49).OdorConcentration = '[10^-7,10^-6]' ;
Input(49).OdorRsp{1} = {'[37:2:50]','[7:2:20]'} ; % no background
Input(49).OdorRsp{2} = {'[38:2:50]','[8:2:20]'} ; % 10^-7 at 100ml/min
Input(49).BgConcentration = {'0','Log7'} ;

% vm7
% +/- background at several pulse concentrations
% small spikes 
Input(50).cellname = '121213_1' ;
Input(50).OdorConcentration = '[10^-7,10^-6,10^-5,10^-4,10^-3]' ;
Input(50).OdorRsp{1} = {'[10:2:23]','[25:2:38]','[41:2:54]','[57:2:70]','[73:2:86]'} ; % no background
Input(50).OdorRsp{2} = {'[11:2:23]','[26:2:38]','[42:2:54]','[58:2:70]','[74:2:86]'} ; % 10^-6 at 100ml/min
Input(50).BgConcentration = {'0','Log6'} ;

% vm7
% +/- background at several pulse concentrations
% small spikes, several technical issues causing delay, several trials show firing rate oscillations
Input(51).cellname = '121220_1' ;
Input(51).OdorConcentration = '[10^-8,10^-7,10^-6,10^-5,10^-4,10^-3]' ;
Input(51).OdorRsp{1} = {'[10:2:23]','[26:2:39]','[41:2:54]','[69:2:78,81:2:88]','[91:2:104]','[107:2:120]'} ; % no background
Input(51).OdorRsp{2} = {'[11:2:23]','[27:2:39]','[42:2:54]','[70:2:78,82:2:88]','[92:2:104]','[108:2:120]'} ; % 10^-5 at 100ml/min
%Input(51).OdorRsp{3} = {'[]','[]','[127:2:140]','[]','[]','[]'} ; % 10^-5 at 100ml/min (did not include control epochs (93:2:103) above)
%Input(51).BgConcentration = {'0','Log5','wash'} ;
Input(51).BgConcentration = {'0','Log5'} ;

% vm7
% +/- background at several pulse concentrations
% poor cell but worth a look (not really usable)
Input(52).cellname = '121227_1' ;
Input(52).OdorConcentration = '[10^-8,10^-7,10^-6,10^-5,10^-4]' ;
Input(52).OdorRsp{1} = {'[6:2:16,20:2:23]','[27:2:30,32:2:41]','[42:2:45]','[46:2:49]','[50:2:53]'} ; % no background
Input(52).OdorRsp{2} = {'[7:2:16,21:2:23]','[28:2:30,33:2:41]','[43:2:45]','[47:2:49]','[51:2:53]'} ; % 10^-5 at 100ml/min
Input(52).BgConcentration = {'0','Log5'} ;

% vm7
% +/- background at several pulse concentrations
% okay cell, small spikes, fly was p4 (all other have been pe2-3)
Input(53).cellname = '121228_1' ;
Input(53).OdorConcentration = '[10^-8,10^-7,10^-6,10^-5,10^-4,10^-3]' ;
Input(53).OdorRsp{1} = {'[4:2:17]','[20:2:33]','[35:2:48]','[51:2:64]','[67:2:80]','[83:2:96]'} ; % no background
Input(53).OdorRsp{2} = {'[5:2:17]','[21:2:33]','[36:2:48]','[52:2:64]','[68:2:80]','[84:2:96]'} ; % 10^-5 at 100ml/min
%Input(53).OdorRsp{3} = {'[]','[]','[103:2:114]','[]','[]','[]'} ; % 10^-5 at 100ml/min (did not include control epochs (102:2:114) above)
% Input(53).BgConcentration = {'0','Log5','wash'} ;
Input(53).BgConcentration = {'0','Log5'} ;

% vm7
% +/- background at several pulse concentrations
% good cell, short data blocks
Input(54).cellname = '121231_1' ;
Input(54).OdorConcentration = '[10^-8,10^-7,10^-6,10^-5,10^-4,10^-3]' ;
Input(54).OdorRsp{1} = {'[4:2:9]','[12:2:17]','[19:2:24]','[27:2:32]','[35:2:40]','[43:2:48]'} ; % no background
Input(54).OdorRsp{2} = {'[5:2:9]','[13:2:17]','[20:2:24]','[28:2:32]','[36:2:40]','[44:2:48]'} ; % 10^-5 at 100ml/min
%Input(54).OdorRsp{3} = {'[]','[]','[52:2:56]','[]','[]','[]'} ; % 10^-5 at 100ml/min (did not include control epochs (51:2:56) above)
% Input(54).BgConcentration = {'0','Log5','wash'} ;
Input(54).BgConcentration = {'0','Log5'} ;

% vm7
% +/- background at several pulse concentrations
% okay cell, small spikes, short data blocks
Input(55).cellname = '130108_1' ;
Input(55).OdorConcentration = '[10^-8,10^-7,10^-6,10^-5,10^-4,10^-3]' ;
Input(55).OdorRsp{1} = {'[4:2:9]','[12:2:17]','[19:2:24]','[27:2:32]','[35:2:40]','[43:2:48]'} ; % no background
Input(55).OdorRsp{2} = {'[5:2:9]','[13:2:17]','[20:2:24]','[28:2:32]','[36:2:40]','[44:2:48]'} ; % 10^-5 at 100ml/min
%Input(55).OdorRsp{3} = {'[]','[]','[52:2:56]','[]','[]','[]'} ; % 10^-5 at 100ml/min (did not include control epochs (51:2:56) above)
%Input(55).BgConcentration = {'0','Log5','wash'} ;
Input(55).BgConcentration = {'0','Log5'} ;

% vm7
% +/- background at several pulse concentrations
% okay cell, small spikes, cell ran down, short data blocks
Input(56).cellname = '130124_1' ;
Input(56).OdorConcentration = '[10^-8,10^-7,10^-6,10^-5]' ; % (data seemingly past run down :'10^-4','10^-3')
Input(56).OdorRsp{1} = {'[6:2:11]','[14:2:19]','[21:2:26]','[29:2:34]'} ; % no background (data seemingly past run down :'[37:2:46]','[49:2:56]')
Input(56).OdorRsp{2} = {'[7:2:11]','[15:2:19]','[22:2:26]','[30:2:34]'} ; % 10^-4 at 100ml/min (data seemingly past run down :'[38:2:46]','[50:2:56]')
Input(56).BgConcentration = {'0','Log4'} ;

% vm7
% +/- background at several pulse concentrations
% good cell, very stable backgrounds, short data blocks, trials at 10^-8 are likely contaminated (response too big)
Input(57).cellname = '130204_1' ;
Input(57).OdorConcentration = '[10^-7,10^-6,10^-5,10^-4,10^-3]' ; %
Input(57).OdorRsp{1} = {'[16:2:23]','[24:2:31]','[34:2:41]','[44:2:51]','[54:2:61]'} ; % no background
Input(57).OdorRsp{2} = {'[17:2:23]','[25:2:31]','[35:2:41]','[45:2:51]','[55:2:61]'} ; % 10^-4 at 100ml/min
Input(57).BgConcentration = {'0','Log4'} ;

% vm7
% +/- background at several pulse concentrations
% okay cell, some variability in rest (but attempted to control with current), short data blocks,  trials at 10^-8 are likely contaminated (response too big)
Input(58).cellname = '130208_1' ;
Input(58).OdorConcentration = '[10^-7,10^-6,10^-5,10^-4,10^-3]' ; %
Input(58).OdorRsp{1} = {'[18:2:25]','[26:2:33]','[36:2:43]','[46:2:53]','[56:2:63]'} ; % no background (cell strongly depolarized after 22)
Input(58).OdorRsp{2} = {'[19:2:25]','[27:2:33]','[37:2:43]','[47:2:53]','[57:2:63]'} ; % 10^-4 at 100ml/min
Input(58).BgConcentration = {'0','Log4'} ;

%  vm7
% +/- background at several pulse concentrations
% not good cell, high variability in rest, short data blocks, trials at 10^-8 are likely contaminated (response too big)
Input(59).cellname = '130213_1' ;
Input(59).OdorConcentration = '[10^-7,10^-6]' ; %
Input(59).OdorRsp{1} = {'[34:2:41]','[45:2:52]'} ; % no background (CAUTION: some trials will need to be excluded)
Input(59).OdorRsp{2} = {'35:2:41]','[46:2:52]'} ; % 10^-7 at 100ml/min
Input(59).BgConcentration = {'0','Log7'} ;

% vm7
% +/- background at several pulse concentrations
% okay cell, some long blocks, used new odor vial caps
Input(60).cellname = '130214_1' ;
Input(60).OdorConcentration = '[10^-8,10^-7,10^-6,10^-5,10^-4,10^-3]' ; %
Input(60).OdorRsp{1} = {'[3:2:10,62:2:71]','[16:2:21,73:2:86]','[22:2:29]','[32:2:39]','[42:2:49]','[52:2:59]'} ; % no background 
Input(60).OdorRsp{2} = {'[4:2:10,63:2:71]','[17:2:21,74:2:86]','[23:2:29]','[33:2:39]','[43:2:49]','[53:2:59]'} ; % 10^-7 at 100ml/min
Input(60).BgConcentration = {'0','Log7'} ;

% vm7
% +/- background at several pulse concentrations
% okay cell, some long blocks
Input(61).cellname = '130218_1' ;
Input(61).OdorConcentration = '[10^-8,10^-7,10^-6,10^-5,10^-4,10^-3]' ; %
Input(61).OdorRsp{1} = {'[4:2:17]','[21:2:34]','[35:2:42]','[45:2:50]','[53:2:58]','[61:2:66]'} ; % no background 
Input(61).OdorRsp{2} = {'[5:2:17]','[22:2:34]','[36:2:42]','[46:2:50]','[54:2:58]','[62:2:66]'} ; % 10^-7 at 100ml/min
Input(61).BgConcentration = {'0','Log7'} ;

% vm7
% +/- background at several pulse concentrations
% okay cell (but lost spikes early), long blocks, more data at 10^-6,-5,and -4 but spike presence is too variable
% CAUTION NORMALIZED DOSE RESPONSE CURVE CANNOT BE TRUSTED (no spike data at hihg odor levels)
Input(62).cellname = '130222_1' ;
Input(62).OdorConcentration = '[10^-8,10^-7]' ; % 
Input(62).OdorRsp{1} = {'[15:2:28]','[31:2:44]'} ; % no background 
Input(62).OdorRsp{2} = {'[16:2:28]','[32:2:44]'} ; % 10^-7 at 100ml/min
Input(62).BgConcentration = {'0','Log7'} ;

% pb1
% +/- background at 2 pulse concentration 
% had very high spontaneous rate and 10^-5 response is very poor
Input(63).cellname = '130424_1' ;
Input(63).OdorConcentration = '[10^-5,10^-4]' ; % 
Input(63).OdorRsp{1} = {'[20:2:24]','[6:2:19]'} ; % no background 
Input(63).OdorRsp{2} = {'[21:2:24]','[7:2:19]'} ; % 10^-5 at 100ml/min
Input(63).BgConcentration = {'0','Log5'} ;

% pb1
% +/- background at 2 pulse concentration 
% had very high spontaneous rate, not sure that this is a different cell than 63 above
Input(64).cellname = '130424_1' ;
Input(64).OdorConcentration = '[10^-5,10^-4]' ; % 
Input(64).OdorRsp{1} = {'[32:2:41]','[42:2:46]'} ; % no background 
Input(64).OdorRsp{2} = {'[33:2:41]','[43:2:46]'} ; % 10^-5 at 100ml/min
Input(64).BgConcentration = {'0','Log5'} ;

% pb1
% +/- background at 1 pulse concentration 
% these spike shapes look different than previous cells(63,64)
Input(65).cellname = '130424_1' ;
Input(65).OdorConcentration = '[10^-5]' ; % 
Input(65).OdorRsp{1} = {'[49:2:55]'} ; % no background 
Input(65).OdorRsp{2} = {'[50:2:55]'} ; % 10^-5 at 100ml/min
Input(65).BgConcentration = {'0','Log5'} ;

% pb1
% +/- background at several pulse concentration 
% spont rate looks high 
% odors a little old and used several times previously: made 4/22
Input(66).cellname = '130426_1' ;
Input(66).OdorConcentration = '[10^-5,10^-4,10^-3]' ; % 
Input(66).OdorRsp{1} = {'[8:2:17]','[18:2:31]','[34:2:47]'} ; % no background 
Input(66).OdorRsp{2} = {'[9:2:17]','[19:2:31]','[35:2:47]'} ; % 10^-5 at 100ml/min
Input(66).BgConcentration = {'0','Log5'} ;

% pb1
% +/- background at 1 pulse concentration 
% cell looked good (had good 10^-5 response, lost when switching odor vial)
Input(67).cellname = '130426_1' ;
Input(67).OdorConcentration = '[10^-7]' ; % 
Input(67).OdorRsp{1} = {'[53:2:82]'} ; % no background 
Input(67).OdorRsp{2} = {'[54:2:82]'} ; % 10^-5 at 100ml/min
Input(67).BgConcentration = {'0','Log5'} ;

% pb1
% +/- background at 1 pulse concentration 
% (lost and regained in this data block)
Input(68).cellname = '130510_1' ;
Input(68).OdorConcentration = '[10^-5]' ; % 
Input(68).OdorRsp{1} = {'[26:2:29,33:2:42]'} ; % no background 
Input(68).OdorRsp{2} = {'[27:2:29,34:2:42]'} ; % 10^-5 at 100ml/min
Input(68).BgConcentration = {'0','Log5'} ;

% pb1
% +/- background at multiple pulse concentration 
% cell looked good (lost and regained a couple times in this data block)
Input(69).cellname = '130510_1' ;
Input(69).OdorConcentration = '[10^-5,10^-4,10^-3]' ; % 
Input(69).OdorRsp{1} = {'[49:2:62]','[63:2:76]','[77:2:84]'} ; % no background 
Input(69).OdorRsp{2} = {'[50:2:62]','[64:2:76]','[78:2:84]'} ; % 10^-5 at 100ml/min
Input(69).BgConcentration = {'0','Log5'} ;

% pb1
% +/- background at single pulse concentration 
% okay cell, spike rate and waveform looked to change 
% epoch 37 lost temp access, control pulse is not visible in psth
Input(70).cellname = '130516_1' ;
Input(70).OdorConcentration = '[10^-6]' ; % 
Input(70).OdorRsp{1} = {'[11:2:46]'} ; % no background 
Input(70).OdorRsp{2} = {'[12:2:46]'} ; % 10^-5 at 100ml/min
Input(70).BgConcentration = {'0','Log5'} ;

% pb1
% +/- background at single pulse concentration 
% sensetive cell but lost access several times 
% did not clear line from previous 10^-5 pulse
Input(71).cellname = '130522_1' ;
Input(71).OdorConcentration = '[10^-6]' ; % 
Input(71).OdorRsp{1} = {'[40:2:45]'} ; % no background 
Input(71).OdorRsp{2} = {'[41:2:45]'} ; % 10^-5 at 100ml/min
Input(71).BgConcentration = {'0','Log5'} ;

% pb1
% +/- background at several pulse concentration 
% good cell, long delay before 10^-7 and odors out of ordor
Input(72).cellname = '130524_1' ;
Input(72).OdorConcentration = '[10^-7,10^-6,10^-5,10^-4,10^-3]' ; % 
Input(72).OdorRsp{1} = {'[107:2:146]','[24:2:29,31:2:52]','[86:2:101]','[55,58:2:67]','[70:2:81]'} ; % no background 
Input(72).OdorRsp{2} = {'[108:2:146]','[25:2:29,32:2:52]','[87:2:101]','[56,59:2:67]','[71:2:81]'} ; % 10^-5 at 100ml/min
Input(72).BgConcentration = {'0','Log5'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at 2 pulse concentration 
% slow spike waveforms
Input(73).cellname = '130530_1' ;
Input(73).OdorConcentration = '[10^-5,10^-3]' ; % 
Input(73).OdorRsp{1} = {'[25:2:34]','[35:2:44]'} ; % no background 
Input(73).OdorRsp{2} = {'[26:2:34]','[36:2:44]'} ; % 10^-5 at 100ml/min
Input(73).BgConcentration = {'0','Log5'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% fly is 5-6 days old (normally 2-3 days old)
% +/- background at several concentration 
% high spont rate during initial 2-but and pCresol pulses but went down during background trials
Input(74).cellname = '130604_1' ;
Input(74).ButPulseRsp = '[19:21]' ; % 2-butanone 10^-5
Input(74).PcPulseRsp = '[22:28]' ; %pCresol 10^-2 
Input(74).OdorConcentration = '[10^-5,10^-3,10^-2]' ; % 
Input(74).OdorRsp{1} = {'[32:2:41]','[44:2:51]','[52:2:57]'} ; % no background 
Input(74).OdorRsp{2} = {'[33:2:41]','[45:2:51]','[53:2:57]'} ; % 10^-5 at 100ml/min
Input(74).BgConcentration = {'0','Log5'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
% lost cell several times with changes in spont rate
Input(75).cellname = '130605_1' ;
Input(75).ButPulseRsp = '[8,10,11,13,14]' ; % 2-butanone 10^-5
Input(75).PcPulseRsp = '[18:24]' ; % pCresol 10^-2 
Input(75).OdorConcentration = '[10^-6,10^-5,10^-4,10^-3,10^-2]' ; % 
Input(75).OdorRsp{1} = {'[31:2:42]','[43:2:48,53:2:60]','[84:2:89,92]','[62:2:69]','[70:2:77]'} ; % no background 
Input(75).OdorRsp{2} = {'[32:2:42]','[44:2:48,54:2:60]','[85:2:89,93]','[63:2:69]','[71:2:77]'} ; % 10^-5 at 100ml/min
Input(75).BgConcentration = {'0','Log5'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
% lost cell several times with changes in spike waveform
Input(76).cellname = '130606_1' ;
Input(76).ButPulseRsp = '[50:2:63]' ; % 2-butanone 10^-6
Input(76).PcPulseRsp = '[64:69]' ; % pCresol 10^-2 
Input(76).OdorConcentration = '[10^-7,10^-6,10^-5,10^-3,10^-2]' ; % 
Input(76).OdorRsp{1} = {'[76:2:93]','[50:2:63]','[15:2:26]','[27:2:34]','[37:2:44]'} ; % no background 
Input(76).OdorRsp{2} = {'[77:2:93]','[51:2:63]','[16:2:26]','[28:2:34]','[38:2:44]'} ; % 10^-5 at 100ml/min
Input(76).BgConcentration = {'0','Log5'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
% lost cell several times 
Input(77).cellname = '130614_1' ;
Input(77).OdorConcentration = '[10^-6,10^-5,10^-4,10^-3,10^-2]' ; % 
Input(77).OdorRsp{1} = {'[72:2:91]','[16:2:31,92:2:95]','[34:2:41]','[44:2:51]','[54,58:2:63]'} ; % no background 
Input(77).OdorRsp{2} = {'[73:2:91]','[17:2:31,93:2:95]','[35:2:41]','[45:2:51]','[55,59:2:63]'} ; % 10^-6 at 100ml/min
Input(77).BgConcentration = {'0','Log6'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
% lost cell several times 
Input(78).cellname = '130621_1' ;
Input(78).OdorConcentration = '[10^-6,10^-5,10^-4,10^-3,10^-2,10^-1]' ; % 
Input(78).OdorRsp{1} = {'[81:2:110]','[27:2:40]','[41:2:48]','[51:2:58]','[61:2:68]','[71:2:76]'} ; % no background 
Input(78).OdorRsp{2} = {'[82:2:110]','[28:2:40]','[42:2:48]','[52:2:58]','[62:2:68]','[72:2:76]'} ; % 10^-6 at 100ml/min
Input(78).BgConcentration = {'0','Log6'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at one pulse concentration 
% lost cell several times 
Input(79).cellname = '130711_1' ;
Input(79).OdorConcentration = '[10^-5]' ; % 
Input(79).OdorRsp{1} = {'[2:2:7]'} ; % no background 
Input(79).OdorRsp{2} = {'[3:2:7]'} ; % 10^-6 at 100ml/min
Input(79).BgConcentration = {'0','Log6'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
% lost cell several times 
Input(80).cellname = '130712_1' ;
Input(80).OdorConcentration = '[10^-5,10^-4,10^-3,10^-2]' ; % 
Input(80).OdorRsp{1} = {'[3,6:2:11,13:2:16]','[19:2:22,24:2:27,28]','[32,35:2:44]','[47:2:56]'} ; % no background 
Input(80).OdorRsp{2} = {'[4,7:2:11,14:2:16]','[20:2:22,25:2:27,29]','[33,36:2:44]','[48:2:56]'} ; % 10^-6 at 100ml/min
Input(80).BgConcentration = {'0','Log6'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
% lost cell several times - 10^-6 and -7 data not included
Input(81).cellname = '130712_1' ;
Input(81).OdorConcentration = '[10^-5,10^-4,10^-3,10^-2]' ; % 
Input(81).OdorRsp{1} = {'[67,72:2:81]','[82:2:91]','[92,95:2:102]','[103:2:112]'} ; % no background 
Input(81).OdorRsp{2} = {'[68,73:2:81]','[83:2:91]','[93,96:2:102]','[104:2:112]'} ; % 10^-6 at 100ml/min
Input(81).BgConcentration = {'0','Log6'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
% lost cell several times 
Input(82).cellname = '130715_1' ;
Input(82).OdorConcentration = '[10^-7,10^-6,10^-5]' ; % 
Input(82).OdorRsp{1} = {'[44:2:63]','[68:2:87]','[21:2:38,88:2:91,93]'} ; % no background 
Input(82).OdorRsp{2} = {'[45:2:63]','[69:2:87]','[22:2:38,89:2:91,94]'} ; % 10^-6 at 100ml/min
Input(82).BgConcentration = {'0','Log6'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
% lost cell several times, very sensitive cell
% some more trials for 10^-7 but spikes detection is poor
Input(83).cellname = '130722_1' ;
Input(83).OdorConcentration = '[10^-7,10^-6,10^-5]' ; % 
Input(83).OdorRsp{1} = {'[59,63]','[35:2:54]','[18:2:23,26]'} ; % no background 
Input(83).OdorRsp{2} = {'[60,64]','[36:2:54]','[19:2:23,27]'} ; % 10^-6 at 100ml/min
Input(83).BgConcentration = {'0','Log6'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
% this ORN seemed to have 2 spike waveforms, second waveform seemed decreased in frequency
Input(84).cellname = '130902_1' ;
Input(84).OdorConcentration = '[10^-5,10^-4,10^-3]' ; % 
Input(84).OdorRsp{1} = {'[17:2:20,23:2:42]','[45:2:52,55:2:64]','[67:2:86]'} ; % no background 
Input(84).OdorRsp{2} = {'[18:2:20,24:2:42]','[46:2:52,56:2:64]','[68:2:86]'} ; % 10^-7 at 100ml/min
Input(84).BgConcentration = {'0','Log7'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at 1 pulse concentration 
% this ORN seemed to have 2 spike waveforms
% more data but spike waveforms and stats seemed to change
Input(85).cellname = '130912_1' ;
Input(85).OdorConcentration = '[10^-5]' ; % 
Input(85).OdorRsp{1} = {'[7:2:12,14]'} ; % no background 
Input(85).OdorRsp{2} = {'[8:2:12,15]'} ; % 10^-7 at 100ml/min
Input(85).BgConcentration = {'0','Log7'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
Input(86).cellname = '130913_1' ;
Input(86).OdorConcentration = '[10^-6,10^-5,10^-4,10^-3]' ; % 
Input(86).OdorRsp{1} = {'[65:2:84,86:2:107]','[47:2:60]','[110:2:129]','[132:2:151]'} ; % no background 
Input(86).OdorRsp{2} = {'[66:2:84,87:2:107]','[48:2:60]','[111:2:129]','[133:2:151]'} ; % 10^-7 at 100ml/min
Input(86).BgConcentration = {'0','Log7'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
Input(87).cellname = '130916_1' ;
Input(87).OdorConcentration = '[10^-6,10^-5,10^-4,10^-3]' ; % 
Input(87).OdorRsp{1} = {'[100:2:113,115:2:146]','[41:2:48,50:2:57,59:2:62]','[63:2:72,74:2:83]','[84:2:97]'} ; % no background 
Input(87).OdorRsp{2} = {'[101:2:113,116:2:146]','[42:2:48,51:2:57,60:2:62]','[64:2:72,75:2:83]','[85:2:97]'} ; % 10^-7 at 100ml/min
Input(87).BgConcentration = {'0','Log7'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
Input(88).cellname = '130918_1' ;
Input(88).OdorConcentration = '[10^-6,10^-5,10^-4,10^-3,10^-2]' ; % 
Input(88).OdorRsp{1} = {'[105:2:150]','[7,10,14:2:29,31:2:34]','[35,41:2:56]','[59:2:70,73:2:78]','[81:2:100]'} ; % no background 
Input(88).OdorRsp{2} = {'[106:2:150]','[8,11,15:2:29,32:2:34]','[36,42:2:56]','[60:2:70,74:2:78]','[82:2:100]'} ; % 10^-4 at 100ml/min
Input(88).BgConcentration = {'0','Log4'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
% WEIRD afterspikes contaminate this data, not sure its a good one to use!
Input(89).cellname = '130919_1' ;
Input(89).OdorConcentration = '[10^-5,10^-4,10^-3,10^-2]' ; % 
Input(89).OdorRsp{1} = {'[26:2:31,33:2:40]','[41:2:48,51:2:62]','[63:2:76]','[77:2:90]'} ; % no background 
Input(89).OdorRsp{2} = {'[27:2:31,34:2:40]','[42:2:48,52:2:62]','[64:2:76]','[78:2:90]'} ; % 10^-4 at 100ml/min
Input(89).BgConcentration = {'0','Log4'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
% second spike contamination looks to be odor responsive (goes away after long pulse)
Input(90).cellname = '130925_1' ;
Input(90).OdorConcentration = '[10^-6,10^-5,10^-4,10^-3,10^-2]' ; % 
Input(90).OdorRsp{1} = {'[101:2:124]','[13:2:32]','[35:2:54]','[57:2:68,70,74:2:77]','[80,84:2:93]'} ; % no background 
Input(90).OdorRsp{2} = {'[102:2:124]','[14:2:32]','[36:2:54]','[58:2:68,71,75:2:77]','[81,85:2:93]'} ; % 10^-4 at 100ml/min
Input(90).BgConcentration = {'0','Log4'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
% difficult detecting spikes during background and many false positives during no odor
Input(91).cellname = '130926_1' ;
Input(91).OdorConcentration = '[10^-6,10^-5,10^-4,10^-3,10^-2]' ; % 
Input(91).OdorRsp{1} = {'[47:2:50,51,54:2:61,63:2:74,77:2:82,87,91]','[25:2:40,43]','[93:2:112]','[113:2:116,119:2:132]','[133:2:152]'} ; % no background 
Input(91).OdorRsp{2} = {'[48:2:50,52,55:2:61,64:2:74,78:2:82,88,92]','[26:2:40,44]','[94:2:112]','[114:2:116,120:2:132]','[134:2:152]'} ; % 10^-4 at 100ml/min
Input(91).BgConcentration = {'0','Log4'} ;

% vm7 (PN - rest of cell in 48)
% +/- background at several pulse concentrations
% small spikes
Input(92).cellname = '121206_1' 
Input(92).OdorConcentration = '[10^-6,10^-5]' ;
Input(92).OdorRsp{1} = {'[7:2:20]','[49:2:62]'} ; % no background
Input(92).OdorRsp{2} = {'[8:2:20]','[50:2:62]'} ; % 10^-7 at 100ml/min
Input(92).BgConcentration = {'0','Log7'} ;

% vm7 (PN- rest of cell 49)
% +/- background at several pulse concentrations
% small spikes and high rest potential
Input(93).cellname = '121212_1' ;
Input(93).OdorConcentration = '[10^-7,10^-6,10^-5,10^-4,10^-3]' ;
Input(93).OdorRsp{1} = {'[51:2:64]','[21:2:34]','[65:2:78]','[79:2:92]','[93:2:106]'} ; % no background
Input(93).OdorRsp{2} = {'[52:2:64]','[22:2:34]','[66:2:78]','[80:2:92]','[94:2:106]'} ; % 10^-6 at 100ml/min
Input(93).BgConcentration = {'0','Log6'} ;

% vm7 (PN- rest of cell 50)
% +/- background at several pulse concentrations
% small spikes 
Input(94).cellname = '121213_1' ;
Input(94).OdorConcentration = '[10^-6]' ;
Input(94).OdorRsp{1} = {'[93:2:102]'} ; % no background
Input(94).OdorRsp{2} = {'[94:2:102]'} ; % 10^-5 at 100ml/min 
Input(94).BgConcentration = {'0','Log5'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
Input(95).cellname = '131004_1' ;
Input(95).OdorConcentration = '[10^-5,10^-4,10^-3]' ; % 
Input(95).OdorRsp{1} = {'[29:2:38]','[41:2:50]','[53:2:62]'} ; % no background 
Input(95).OdorRsp{2} = {'[30:2:38]','[42:2:50]','[54:2:62]'} ; % 10^-7 at 100ml/min
Input(95).BgConcentration = {'0','Log7'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
Input(96).cellname = '131004_1' ;
Input(96).OdorConcentration = '[10^-6,10^-5]' ; % 
Input(96).OdorRsp{1} = {'[91:2:126]','[67:2:86]'} ; % no background 
Input(96).OdorRsp{2} = {'[92:2:126]','[68:2:86]'} ; % 10^-7 at 100ml/min
Input(96).BgConcentration = {'0','Log7'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
Input(97).cellname = '131016_1' ;
Input(97).OdorConcentration = '[10^-6,10^-5,10^-4]' ; % 
Input(97).OdorRsp{1} = {'[30:2:69]','[15:2:24]','[71:2:92]'} ; % no background 
Input(97).OdorRsp{2} = {'[31:2:69]','[16:2:24]','[72:2:92]'} ; % 10^-7 at 100ml/min
Input(97).BgConcentration = {'0','Log7'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
Input(98).cellname = '131017_1' ;
Input(98).OdorConcentration = '[10^-6,10^-3,10^-2]' ; % 
Input(98).OdorRsp{1} = {'[13:2:22,24:2:49]','[50:2:55,58:2:69]','[70:2:89]'} ; % no background 
Input(98).OdorRsp{2} = {'[14:2:22,25:2:49]','[51:2:55,59:2:69]','[71:2:89]'} ; % 10^-7 at 100ml/min
Input(98).BgConcentration = {'0','Log7'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
Input(99).cellname = '131018_1' ;
Input(99).OdorConcentration = '[10^-6,10^-5,10^-4,10^-3,10^-2]' ; % 
Input(99).OdorRsp{1} = {'[62:2:101]','[12,15:2:24]','[27:2:40]','[43:2:56]','[104:2:111]'} ; % no background 
Input(99).OdorRsp{2} = {'[63:2:101]','[13,16:2:24]','[28:2:40]','[44:2:56]','[105:2:111]'} ; % 10^-6 at 100ml/min
Input(99).BgConcentration = {'0','Log6'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
% some useable trials at 10^-2 but not many and hard to get propper spike data
Input(100).cellname = '131030_1' ;
Input(100).OdorConcentration = '[10^-6,10^-5,10^-4]' ; % 
Input(100).OdorRsp{1} = {'[32:2:37,40:2:69]','[20:2:27]','[72:2:79,82]'} ; % no background 
Input(100).OdorRsp{2} = {'[33:2:37,41:2:69]','[21:2:27]','[73:2:79,83]'} ; % 10^-5 at 100ml/min
Input(100).BgConcentration = {'0','Log5'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
Input(101).cellname = '131031_1' ;
Input(101).OdorConcentration = '[10^-5,10^-2]' ; % 
Input(101).OdorRsp{1} = {'[27:2:32]','[33:2:42,44:2:51]'} ; % no background 
Input(101).OdorRsp{2} = {'[28:2:32]','[34:2:42,45:2:51]'} ; % 10^-5 at 100ml/min
Input(101).BgConcentration = {'0','Log5'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
Input(102).cellname = '131031_1' ;
Input(102).OdorConcentration = '[10^-2]' ; % 
Input(102).OdorRsp{1} = {'[53:2:72]'} ; % no background 
Input(102).OdorRsp{2} = {'[54:2:72]'} ; % 10^-5 at 100ml/min
Input(102).BgConcentration = {'0','Log5'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
% data available at 10^-5([8:21]) but spike detection needs improvment to extract spike times
Input(103).cellname = '131101_1' ;
Input(103).OdorConcentration = '[10^-6,10^-4,10^-3,10^-2]' ; % 
Input(103).OdorRsp{1} = {'[27:2:42,44:2:63]','[66:2:73,76:2:85]','[88:2:103]','[106:2:125]'} ; % no background 
Input(103).OdorRsp{2} = {'[28:2:42,45:2:63]','[67:2:73,77:2:85]','[89:2:103]','[107:2:125]'} ; % 10^-4 at 100ml/min
Input(103).BgConcentration = {'0','Log4'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
% Spike detection is NOT GOOD in this cell (Might need to rework code to avoid these sporadic issues)
Input(104).cellname = '131106_1' ;
Input(104).OdorConcentration = '[10^-6,10^-5,10^-4,10^-3,10^-2]' ; % 
Input(104).OdorRsp{1} = {'[28,31:2:34,37:2:46,48:2:53,55:2:74]','[12:2:23]','[75:2:94]','[95:2:114]','[115:2:124]'} ; % no background 
Input(104).OdorRsp{2} = {'[29,32:2:34,38:2:46,49:2:53,56:2:74]','[13:2:23]','[76:2:94]','[96:2:114]','[116:2:124]'} ; % 10^-4 at 100ml/min
Input(104).BgConcentration = {'0','Log4'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
Input(105).cellname = '131107_1' ;
Input(105).OdorConcentration = '[10^-6,10^-5]' ; % 
Input(105).OdorRsp{1} = {'[30:2:69]','[15:2:20,22:2:25]'} ; % no background 
Input(105).OdorRsp{2} = {'[31:2:69]','[16:2:20,23:2:25]'} ; % 10^-6 at 100ml/min
Input(105).BgConcentration = {'0','Log6'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
% same cell as 105 
% trial 87 has a anlomalous post pulse response, but is left becease it does not impact any important stats
Input(106).cellname = '131107_1' ;
Input(106).OdorConcentration = '[10^-2]' ; % 
Input(106).OdorRsp{1} = {'[70:2:89]'} ; % no background 
Input(106).OdorRsp{2} = {'[71:2:89]'} ; % 10^-5 at 100ml/min
Input(106).BgConcentration = {'0','Log5'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
Input(107).cellname = '131108_1' ;
Input(107).OdorConcentration = '[10^-2]' ; % 
Input(107).OdorRsp{1} = {'[10:2:13,16:2:19,21:2:34]'} ; % no background 
Input(107).OdorRsp{2} = {'[11:2:13,17:2:19,22:2:34]'} ; % 10^-7 at 100ml/min
Input(107).BgConcentration = {'0','Log7'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
Input(108).cellname = '131108_1' ;
Input(108).OdorConcentration = '[10^-2]' ; % 
Input(108).OdorRsp{1} = {'[43,48:2:51,53:2:66]'} ; % no background 
Input(108).OdorRsp{2} = {'[44,49:2:51,54:2:66]'} ; % 10^-7 at 100ml/min
Input(108).BgConcentration = {'0','Log7'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
Input(109).cellname = '131111_1' ;
Input(109).OdorConcentration = '[10^-2]' ; % 
Input(109).OdorRsp{1} = {'[5,8,10]'} ; % no background 
Input(109).OdorRsp{2} = {'[6,9,11]'} ; % 10^-7 at 100ml/min
Input(109).BgConcentration = {'0','Log7'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
Input(110).cellname = '131111_1' ;
Input(110).OdorConcentration = '[10^-2]' ; % 
Input(110).OdorRsp{1} = {'[33:2:40]'} ; % no background 
Input(110).OdorRsp{2} = {'[34:2:40]'} ; % 10^-7 at 100ml/min
Input(110).BgConcentration = {'0','Log7'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
Input(111).cellname = '131111_1' ;
Input(111).OdorConcentration = '[10^-2]' ; % 
Input(111).OdorRsp{1} = {'[57:2:62,64:2:69]'} ; % no background 
Input(111).OdorRsp{2} = {'[58:2:62,65:2:69]'} ; % 10^-7 at 100ml/min
Input(111).BgConcentration = {'0','Log7'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% +/- background at several concentration 
Input(112).cellname = '131111_1' ;
Input(112).OdorConcentration = '[10^-2]' ; % 
Input(112).OdorRsp{1} = {'[73:2:76,78:2:81,83:2:92]'} ; % no background 
Input(112).OdorRsp{2} = {'[74:2:76,79:2:81,84:2:92]'} ; % 10^-7 at 100ml/min
Input(112).BgConcentration = {'0','Log7'} ;

% *IMPORTANT CHANGE BELOW****************
% data blocks above (46-112) are using 2-butanone 
% below are cells using same paradym using ethyl acetate

% pb1 (fly or71a mutant 0A9 "T2")
% ethyl acetate pulse +/- background
Input(113).cellname = '131115_1' ;
Input(113).OdorConcentration = '[10^-4]' ; % 
Input(113).OdorRsp{1} = {'[3:2:12,14:2:17,19:2:22,24]'} ; % no background 
Input(113).OdorRsp{2} = {'[4:2:12,15:2:17,20:2:22,25]'} ; % 10^-5 at 100ml/min
Input(113).BgConcentration = {'0','Log5'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% ethyl acetate pulse +/- background
% same cell as 113
Input(114).cellname = '131115_1' ;
Input(114).OdorConcentration = '[10^-3]' ; % 
Input(114).OdorRsp{1} = {'[26:2:33,35:2:50]'} ; % no background 
Input(114).OdorRsp{2} = {'[27:2:33,36:2:50]'} ; % 10^-4 at 100ml/min
Input(114).BgConcentration = {'0','Log4'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% ethyl acetate pulse +/- background
Input(115).cellname = '131118_1' ;
Input(115).OdorConcentration = '[10^-4]' ; % 
Input(115).OdorRsp{1} = {'[17:2:22,25:2:34,36:2:41]'} ; % no background 
Input(115).OdorRsp{2} = {'[18:2:22,26:2:34,37:2:41]'} ; % 10^-5 at 100ml/min
Input(115).BgConcentration = {'0','Log5'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% ethyl acetate pulse +/- background
Input(116).cellname = '131121_1' ;
Input(116).OdorConcentration = '[10^-4]' ; % 
Input(116).OdorRsp{1} = {'[40:2:45,48:2:55,58:2:77]'} ; % no background 
Input(116).OdorRsp{2} = {'[41:2:45,49:2:55,59:2:77]'} ; % 10^-5 at 100ml/min
Input(116).BgConcentration = {'0','Log5'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% ethyl acetate pulse +/- background
% same cell as 116
Input(117).cellname = '131121_1' ;
Input(117).OdorConcentration = '[10^-3]' ; % 
Input(117).OdorRsp{1} = {'[78:2:107]'} ; % no background 
Input(117).OdorRsp{2} = {'[79:2:107]'} ; % 10^-4 at 100ml/min
Input(117).BgConcentration = {'0','Log4'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% ethyl acetate pulse +/- background
% background was left with tubes unconnected over weekend (exclude?)
Input(118).cellname = '131125_1' ;
Input(118).OdorConcentration = '[10^-4]' ; % 
Input(118).OdorRsp{1} = {'[11,15,17,20:2:35]'} ; % no background 
Input(118).OdorRsp{2} = {'[12,16,18,21:2:35]'} ; % 10^-5 at 100ml/min
Input(118).BgConcentration = {'0','Log5'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% ethyl acetate pulse +/- background
% same cell as 118
Input(119).cellname = '131125_1' ;
Input(119).OdorConcentration = '[10^-3]' ; % 
Input(119).OdorRsp{1} = {'[36:2:75]'} ; % no background 
Input(119).OdorRsp{2} = {'[37:2:75]'} ; % 10^-4 at 100ml/min
Input(119).BgConcentration = {'0','Log4'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% ethyl acetate pulse +/- background
% more trials but getting too many double spikes (90:119)
% no real pulse response detected in accepted trials 
Input(120).cellname = '131125_1' ;
Input(120).OdorConcentration = '[10^-4]' ; % 
Input(120).OdorRsp{1} = {'[90,92,102,110,114:2:118]'} ; % no background 
Input(120).OdorRsp{2} = {'[91,93,103,111,115:2:119]'} ; % 10^-5 at 100ml/min
Input(120).BgConcentration = {'0','Log5'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% ethyl acetate pulse +/- background
% same cell as 120
Input(121).cellname = '131125_1' ;
Input(121).OdorConcentration = '[10^-3]' ; % 
Input(121).OdorRsp{1} = {'[120:2:139,143:2:150]'} ; % no background 
Input(121).OdorRsp{2} = {'[121:2:139,144:2:150]'} ; % 10^-4 at 100ml/min
Input(121).BgConcentration = {'0','Log4'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% ethyl acetate pulse +/- background
Input(122).cellname = '131127_1' ;
Input(122).OdorConcentration = '[10^-4]' ; % 
Input(122).OdorRsp{1} = {'[33:2:36,39:2:42,45,48]'} ; % no background 
Input(122).OdorRsp{2} = {'[34:2:36,40:2:42,46,49]'} ; % 10^-5 at 100ml/min
Input(122).BgConcentration = {'0','Log5'} ;

% pb1 (fly or71a mutant 0A9 "T2")
% ethyl acetate pulse +/- background
% same cell as 122
Input(123).cellname = '131127_1' ;
Input(123).OdorConcentration = '[10^-3]' ; % 
Input(123).OdorRsp{1} = {'[50:2:53,55,57,60:2:69,72:2:85]'} ; % no background 
Input(123).OdorRsp{2} = {'[51:2:53,56,58,61:2:69,73:2:85]'} ; % 10^-4 at 100ml/min
Input(123).BgConcentration = {'0','Log4'} ;

% *BACK TO PNS BELOW

% vm7 
% +/- background at several pulse concentrations
% antennas wet and missing 
Input(124).cellname = '131216_1' ;
Input(124).OdorConcentration = '[10^-6,10^-5]' ;
Input(124).OdorRsp{1} = {'[3:2:16]','[19:2:32]'} ; % no background
Input(124).OdorRsp{2} = {'[4:2:16]','[20:2:32]'} ; % 10^-6 at 100ml/min 
Input(124).BgConcentration = {'0','Log6'} ;

% vm7 
% +/- background at several pulse concentrations
% very small spikes 
Input(125).cellname = '131216_1' ;
Input(125).OdorConcentration = '[10^-6,10^-4,10^-3]' ;
Input(125).OdorRsp{1} = {'[45:2:54]','[57:2:76]','[79:2:96]'} ; % no background
Input(125).OdorRsp{2} = {'[46:2:54]','[58:2:76]','[80:2:96]'} ; % 10^-6 at 100ml/min 
Input(125).BgConcentration = {'0','Log6'} ;

% vm7 
% +/- background at several pulse concentrations
Input(126).cellname = '131219_1' ;
Input(126).OdorConcentration = '[10^-7,10^-6,10^-5,10^-4,10^-3]' ;
Input(126).OdorRsp{1} = {'[60:2:73]','[3:2:12]','[15:2:24]','[27:2:36]','[39:2:52]'} ; % no background
Input(126).OdorRsp{2} = {'[61:2:73]','[4:2:12]','[16:2:24]','[28:2:36]','[40:2:52]'} ; % 10^-4 at 100ml/min 
Input(126).BgConcentration = {'0','Log4'} ;

% vm7 
% +/- background at several pulse concentrations
Input(127).cellname = '131230_1' ;
Input(127).OdorConcentration = '[10^-7,10^-6,10^-5,10^-4,10^-3]' ;
Input(127).OdorRsp{1} = {'[49:2:62]','[5:2:14]','[15:2:24]','[25:2:34]','[35:2:44]'} ; % no background
Input(127).OdorRsp{2} = {'[50:2:62]','[6:2:14]','[16:2:24]','[26:2:34]','[36:2:44]'} ; % 10^-4 at 100ml/min 
Input(127).BgConcentration = {'0','Log4'} ;

% vm7 
% +/- background at several pulse concentrations
% very small spikes  
Input(128).cellname = '140102_1' ;
Input(128).OdorConcentration = '[10^-6,10^-5]' ;
Input(128).OdorRsp{1} = {'[7:2:16]','[19:2:32]'} ; % no background
Input(128).OdorRsp{2} = {'[8:2:16]','[20:2:32]'} ; % 10^-6 at 100ml/min 
Input(128).BgConcentration = {'0','Log6'} ;

% vm7 
% +/- background at several pulse concentrations
% a lot of variablility - may not be a keeper 
Input(129).cellname = '140103_1' ;
Input(129).OdorConcentration = '[10^-6]' ;
Input(129).OdorRsp{1} = {'[8:2:15]',} ; % no background
Input(129).OdorRsp{2} = {'[9:2:15]',} ; % 10^-6 at 100ml/min 
Input(129).BgConcentration = {'0','Log6'} ;

% vm7 
% +/- background at several pulse concentrations
% spikes before trial 30 may not be detectable
Input(130).cellname = '140109_1' ;
Input(130).OdorConcentration = '[10^-7,10^-4,10^-3]' ;
Input(130).OdorRsp{1} = {'[30:2:43]','[46:2:61]','[64:2:77]'} ; % no background
Input(130).OdorRsp{2} = {'[31:2:43]','[47:2:61]','[65:2:77]'} ; % 10^-6 at 100ml/min 
Input(130).BgConcentration = {'0','Log6'} ;

% vm7 
% +/- background at several pulse concentrations
% spikes may not be detectable after trial 13
Input(131).cellname = '140110_1' ;
Input(131).OdorConcentration = '[10^-7,10^-6,10^-5,10^-4,10^-3]' ;
Input(131).OdorRsp{1} = {'[18:2:31]','[32:2:45]','[48:2:61]','[62:2:75]','[2:2:13]'} ; % no background
Input(131).OdorRsp{2} = {'[19:2:31]','[33:2:45]','[49:2:61]','[63:2:75]','[3:2:13]'} ; % 10^-6 at 100ml/min 
Input(131).BgConcentration = {'0','Log6'} ;

%vm7
% +/- background at 1 pulse concentration
% spikes are not reliably detectable
Input(132).cellname = '140203_1' ;
Input(132).OdorConcentration = '[10^-4]' ;
Input(132).OdorRsp{1} = {'[7:2:14]'} ; % no background
Input(132).OdorRsp{2} = {'[8:2:14]'} ; % 10^-7 at 100ml/min 
Input(132).BgConcentration = {'0','Log7'} ;

% vm7 
% +/- background at several pulse concentrations
% 10^-5 data collected but spike detection is difficult (prob not impossible)
Input(133).cellname = '140207_1' ;
Input(133).OdorConcentration = '[10^-7,10^-6,10^-4,10^-3]' ;
Input(133).OdorRsp{1} = {'[42:2:53]','[55:2:64]','[3:2:16]','[19:2:32]'} ; % no background
Input(133).OdorRsp{2} = {'[43:2:53]','[56:2:64]','[4:2:16]','[20:2:32]'} ; % 10^-7 at 100ml/min 
Input(133).BgConcentration = {'0','Log7'} ;

% vm7 
% +/- background at several pulse concentrations
Input(134).cellname = '140210_1' ;
Input(134).OdorConcentration = '[10^-7,10^-6,10^-5,10^-4,10^-3]' ;
Input(134).OdorRsp{1} = {'[25:2:34]','[35:2:42]','[43:2:56]','[2:2:11]','[12:2:21]'} ; % no background
Input(134).OdorRsp{2} = {'[26:2:34]','[36:2:42]','[44:2:56]','[3:2:11]','[13:2:21]'} ; % 10^-7 at 100ml/min 
Input(134).BgConcentration = {'0','Log7'} ;

% vm7 
% +/- background at several pulse concentrations
% spike detection is a bit too sensitive 
Input(135).cellname = '140214_1' ;
Input(135).OdorConcentration = '[10^-7,10^-6,10^-5,10^-4,10^-3]' ;
Input(135).OdorRsp{1} = {'[64:2:79]','[14:2:21]','[24:2:33]','[36:2:45]','[48:2:57]'} ; % no background
Input(135).OdorRsp{2} = {'[65:2:79]','[15:2:21]','[25:2:33]','[37:2:45]','[49:2:57]'} ; % 10^-5 at 100ml/min 
Input(135).BgConcentration = {'0','Log5'} ;

% vm7 
% +/- background at one pulse concentrations
% cell was not in great shape
Input(136).cellname = '140217_1' ;
Input(136).OdorConcentration = '[10^-6]' ;
Input(136).OdorRsp{1} = {'[2:2:11]'} ; % no background
Input(136).OdorRsp{2} = {'[3:2:11]'} ; % 10^-5 at 100ml/min 
Input(136).BgConcentration = {'0','Log5'} ;

% vm7 
% +/- background at several pulse concentrations
% data at lower pulse but spiking seemed substantially depressed
Input(137).cellname = '140220_1' ;
Input(137).OdorConcentration = '[10^-4,10^-3]' ;
Input(137).OdorRsp{1} = {'[5:2:16]','[19:2:30]'} ; % no background
Input(137).OdorRsp{2} = {'[6:2:16]','[20:2:30]'} ; % 10^-4 at 100ml/min 
Input(137).BgConcentration = {'0','Log4'} ;

% vm7 
% +/- background at several pulse concentrations
Input(138).cellname = '140226_1' ;
Input(138).OdorConcentration = '[10^-7,10^-6,10^-5,10^-4,10^-3]' ;
Input(138).OdorRsp{1} = {'[13:2:20,23:2:28]','[33,37:2:42]','[45:2:54]','[57:2:64]','[71:2:82]'} ; % no background
Input(138).OdorRsp{2} = {'[14:2:20,24:2:28]','[34,38:2:42]','[46:2:54]','[58:2:64]','[72:2:83]'} ; % 10^-6 at 100ml/min 
Input(138).BgConcentration = {'0','Log6'} ;

% vm7 
% +/- background at several pulse concentrations
Input(139).cellname = '140227_1' ;
Input(139).OdorConcentration = '[10^-6,10^-5,10^-4,10^-3]' ;
Input(139).OdorRsp{1} = {'[13:2:24]','[27:2:36]','[39:2:48]','[51:2:64]'} ; % no background
Input(139).OdorRsp{2} = {'[14:2:24]','[28:2:36]','[40:2:48]','[52:2:64]'} ; % 10^-5 at 100ml/min 
Input(139).BgConcentration = {'0','Log5'} ;


% ***ethyl acetate PNs below 

% vm7 
% ethyl acetate pulse +/- background
Input(140).cellname = '140228_1' ;
Input(140).OdorConcentration = '[10^-4]' ;
Input(140).OdorRsp{1} = {'[9:2:24]'} ; % no background
Input(140).OdorRsp{2} = {'[10:2:24]'} ; % 10^-5 at 100ml/min 
Input(140).BgConcentration = {'0','Log5'} ;

% vm7 
% ethyl acetate pulse +/- background
Input(141).cellname = '140228_1' ;
Input(141).OdorConcentration = '[10^-3]' ;
Input(141).OdorRsp{1} = {'[27:2:50]'} ; % no background
Input(141).OdorRsp{2} = {'[28:2:50]'} ; % 10^-4 at 100ml/min 
Input(141).BgConcentration = {'0','Log4'} ;

% vm7 
% cell had frequent occurances of extended (1-3+ secs) hyperpolarization that should be excluded
% ethyl acetate pulse +/- background
Input(142).cellname = '140313_1' ;
Input(142).OdorConcentration = '[10^-4]' ;
Input(142).OdorRsp{1} = {'[10:2:31]'} ; % no background
Input(142).OdorRsp{2} = {'[11:2:31]'} ; % 10^-5 at 100ml/min 
Input(142).BgConcentration = {'0','Log5'} ;

% vm7 
% cell had frequent occurances of extended (1-3+ secs) hyperpolarization that should be excluded
% ethyl acetate pulse +/- background
Input(143).cellname = '140313_1' ;
Input(143).OdorConcentration = '[10^-3]' ;
Input(143).OdorRsp{1} = {'[32:2:51]'} ; % no background
Input(143).OdorRsp{2} = {'[33:2:51]'} ; % 10^-4 at 100ml/min 
Input(143).BgConcentration = {'0','Log4'} ;

% vm7 
% ethyl acetate pulse +/- background
Input(144).cellname = '140317_1' ;
Input(144).OdorConcentration = '[10^-4]' ;
Input(144).OdorRsp{1} = {'[4:2:19]'} ; % no background
Input(144).OdorRsp{2} = {'[5:2:19]'} ; % 10^-5 at 100ml/min 
Input(144).BgConcentration = {'0','Log5'} ;

% vm7 
% ethyl acetate pulse +/- background
Input(145).cellname = '140317_1' ;
Input(145).OdorConcentration = '[10^-3]' ;
Input(145).OdorRsp{1} = {'[20:2:39]'} ; % no background
Input(145).OdorRsp{2} = {'[21:2:39]'} ; % 10^-4 at 100ml/min 
Input(145).BgConcentration = {'0','Log4'} ;

% vm7 
% ethyl acetate pulse +/- background
Input(146).cellname = '140321_1' ;
Input(146).OdorConcentration = '[10^-4]' ;
Input(146).OdorRsp{1} = {'[6:2:25]'} ; % no background
Input(146).OdorRsp{2} = {'[7:2:25]'} ; % 10^-5 at 100ml/min 
Input(146).BgConcentration = {'0','Log5'} ;

% vm7 
% ethyl acetate pulse +/- background
Input(147).cellname = '140321_1' ;
Input(147).OdorConcentration = '[10^-3]' ;
Input(147).OdorRsp{1} = {'[28:2:47]'} ; % no background
Input(147).OdorRsp{2} = {'[29:2:47]'} ; % 10^-4 at 100ml/min 
Input(147).BgConcentration = {'0','Log4'} ;

% vm7 
% ethyl acetate pulse +/- background
Input(148).cellname = '140321_1' ;
Input(148).OdorConcentration = '[10^-4]' ;
Input(148).OdorRsp{1} = {'[58:2:77]'} ; % no background
Input(148).OdorRsp{2} = {'[59:2:77]'} ; % 10^-5 at 100ml/min 
Input(148).BgConcentration = {'0','Log5'} ;

% vm7 
% ethyl acetate pulse +/- background
Input(149).cellname = '140321_1' ;
Input(149).OdorConcentration = '[10^-3]' ;
Input(149).OdorRsp{1} = {'[78:2:97]'} ; % no background
Input(149).OdorRsp{2} = {'[79:2:97]'} ; % 10^-4 at 100ml/min 
Input(149).BgConcentration = {'0','Log4'} ;

%********** changed odor pulse valve 3/27
%********** testing odor with PID below (see notes on odor prep 3/27/14)

%PID (gain=x10, pump=high, offset=9+)
% There are plain parafin oil trials as well and some trials where 2nd valve line for bg was not attached
% 2-butanone +/-background at several pulse concentrations
Input(150).cellname = '140327_1' ;
Input(150).OdorConcentration = '[10^-7,10^-6,10^-5,10^-4,10^-3,10^-2]' ;
Input(150).OdorRsp{1} = {'[72:2:91]','[93:2:112]','[113:2:132]','[135:2:154]','[155:164]','[166:175]'} ; % no background
Input(150).OdorRsp{2} = {'[73:2:91]','[]','[]','[]','[]','[]'} ; % 10^-7 at 100ml/min 
Input(150).OdorRsp{3} = {'[]','[94:2:112]','[]','[]','[]','[]'} ; % 10^-6 at 100ml/min 
Input(150).OdorRsp{4} = {'[]','[]','[114:2:132]','[]','[]','[]'} ; % 10^-5 at 100ml/min 
Input(150).OdorRsp{5} = {'[]','[]','[]','[136:2:154]','[]','[]'} ; % 10^-4 at 100ml/min 
Input(150).BgConcentration = {'0','Log7','Log6','Log5','Log4'} ;

%PID (gain=x10, pump=high, offset=9+)
% 10^-9 trials are plain parafin oil
% There are 10^-8 trials, and no vial trials as well
% there are trials where background pulse bottle and vial were switched
% 2-butanone +/-background at several pulse concentrations
Input(151).cellname = '140424_1' ;
Input(151).OdorConcentration = '[10^-9,10^-7,10^-6,10^-5,10^-4,10^-3,10^-2]' ;
Input(151).OdorRsp{1} = {'[10:15,22:25]','[39:2:58]','[59:2:78]','[79:2:98]','[99:2:118]','[119:2:138]','[139:2:158]'} ; % no background
Input(151).OdorRsp{2} = {'[]','[40:2:58]','[]','[]','[]','[]','[]'} ; % 10^-7 at 100ml/min 
Input(151).OdorRsp{3} = {'[]','[]','[60:2:78]','[]','[]','[]','[]'} ; % 10^-6 at 100ml/min 
Input(151).OdorRsp{4} = {'[]','[]','[]','[80:2:98]','[]','[]','[]'} ; % 10^-5 at 100ml/min 
Input(151).OdorRsp{5} = {'[]','[163:2:181]','[183:2:201]','[203:2:221]','[100:2:118]','[120:2:138]','[140:2:158]'} ; % 10^-4 at 100ml/min 
Input(151).BgConcentration = {'0','Log7','Log6','Log5','Log4'} ;

%PID (gain=x10, pump=high, offset=9+)
% This is the same odor bottles as data block 151
% 2-butanone +/-background at several pulse concentrations
Input(152).cellname = '140424_1' ;
Input(152).OdorConcentration = '[10^-7,10^-6,10^-5,10^-4,10^-3,10^-2]' ;
Input(152).OdorRsp{1} = {'[162:2:181]','[182:2:201]','[202:2:221]','[99:2:118]','[119:2:138]','[139:2:158]'} ; % no background
Input(152).OdorRsp{2} = {'[163:2:181]','[183:2:201]','[203:2:221]','[100:2:118]','[120:2:138]','[140:2:158]'} ; % 10^-4 at 100ml/min 
Input(152).BgConcentration = {'0','Log4'} ;



%% prep for export to hdf5 file and analyze

% prep matrix if it doesn't exist already
if ~exist('ForIgor')
    ForIgor=struct() ;
end

% analysis functions

% analysis of PN odor response +/- background

if perform.PNadaptationAnalysisI ;
    Temp = PNadaptationAnalysisI(Input,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.PNadaptationAnalysisJ ;
    Temp = PNadaptationAnalysisJ(Input,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.PNadaptationAnalysisK ;
    Temp = PNadaptationAnalysisK(Input,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.PNadaptationAnalysisL ;
    Temp = PNadaptationAnalysisL(Input,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.adaptationAnalysisApn ;
    Temp = adaptationAnalysisA(Input,A,'PNdata')
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.adaptationAnalysisBpn ;
    Temp = adaptationAnalysisB(Input,A,'PNdata')
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.adaptationAnalysisDpn ;
    Temp = adaptationAnalysisD(Input,A,'PNdata')
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.adaptationAnalysisEpn ;
    Temp = adaptationAnalysisE(Input,A,'PNdata')
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.adaptationAnalysisFpn ;
    Temp = adaptationAnalysisF(Input,A,'PNdata')
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

% analysis of ORN odor response +/- background
if perform.ORNadaptationAnalysisA ;
    Temp = ORNadaptationAnalysisA(Input,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.ORNadaptationAnalysisB ;
    Temp = ORNadaptationAnalysisB(Input,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.ORNadaptationAnalysisC ;
    Temp = ORNadaptationAnalysisC(Input,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.adaptationAnalysisAorn ;
    Temp = adaptationAnalysisA(Input,A,'ORNdata')
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.adaptationAnalysisBorn ;
    Temp = adaptationAnalysisB(Input,A,'ORNdata')
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.adaptationAnalysisDorn ;
    Temp = adaptationAnalysisD(Input,A,'ORNdata')
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.adaptationAnalysisEorn ;
    Temp = adaptationAnalysisE(Input,A,'ORNdata')
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.adaptationAnalysisForn ;
    Temp = adaptationAnalysisF(Input,A,'ORNdata')
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

% analysis of PID data
if perform.PIDAnalysisA ;
    Temp = PIDAnalysisA(Input,A,'')
    ForIgor = mergeStruct(ForIgor,Temp) ;
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

% if max(fieldLength)<32 ;
%     cd Z:/cafaro' Documents'/Analysis/ForIgorHdfs/   
%     exportStructToHDF5(ForIgor,'PNadaptationAnalysisI.h5','/')
% else disp('name too long for igor')
% end




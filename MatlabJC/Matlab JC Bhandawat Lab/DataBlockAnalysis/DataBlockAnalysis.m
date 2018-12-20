%Data notes, data blocks and analysis logicals
%Bhandawat lab
%J Cafaro 1/4/12

%% Analysis logicals (true or false) and data blocks to be analyzed
perform.ORNanalysisA = false ; % this one is h5 data extraction

perform.ORNanalysisB = false ; % this one is h5 data extraction

perform.PNadaptationAnalysisA = false ; 

perform.PNadaptationAnalysisB = false ; % non-existent

perform.ValveNoiseAnalysis = false ;

perform.PNadaptationAnalysisC = false; % this looks at PN adaptation +/- an odor background during a set of odor pulses 

perform.PNadaptationAnalysisD = false ; % similar to C but better parameter timing consistency and logic

perform.PNadaptationAnalysisE = false ; % similar to D but changed odor rsp time assessment 

perform.PNadaptationAnalysisF = false ; % similar to D but changed odor rsp time assessment 

perform.PNadaptationAnalysisG = false ; % similar to F but set up to look at double valve experiments (i.e. background is on valve opened and closed during each trial)

perform.PNadaptationAnalysisH = true ; % similar to G but using different time points for voltage and spike data parameters

perform.ValveNoiseAnalysisB = false ; % this looks at adaptation during a valve noise stimulus

% perform analysis of following data blocks...
for A = [46:51,53:55] ;
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
Input(46).cellname = '121203_1' ;
Input(46).OdorConcentration = '[10^-7,10^-6,10^-5,10^-4,10^-3]' ;
Input(46).OdorRsp{1} = {'[11:2:24]','[25:2:42]','[60:2:75]','[43:2:56]','[76:2:91]'} ; % no background
Input(46).OdorRsp{2} = {'[12:2:24]','[26:2:42]','[61:2:75]','[44:2:56]','[77:2:91]'} ; % 10^-7 at 100ml/min
Input(46).BgConcentration = {'0','Log7'} ;

% VM7
% +/- background at several pulse concentrations
Input(47).cellname = '121205_1' 
Input(47).OdorConcentration = '[10^-7,10^-6,10^-5,10^-4,10^-3]' ;
Input(47).OdorRsp{1} = {'[4:2:17]','[20:2:33]','[36:2:49]','[52:2:65]','[68:2:73]'} ; % no background
Input(47).OdorRsp{2} = {'[5:2:17]','[21:2:33]','[37:2:49]','[53:2:65]','[69:2:73]'} ; % 10^-7 at 100ml/min
Input(47).BgConcentration = {'0','Log7'} ;

% vm7
% +/- background at several pulse concentrations
% small spikes
Input(48).cellname = '121206_1' 
Input(48).OdorConcentration = '[10^-6,10^-5,10^-4]' ;
Input(48).OdorRsp{1} = {'[7:2:20,21:2:34]','[35:2:48,49:2:62]','[63:2:76]'} ; % no background
Input(48).OdorRsp{2} = {'[8:2:20]','[50:2:62]','[]'} ; % 10^-7 at 100ml/min
Input(48).OdorRsp{3} = {'[22:2:34]','[36:2:48]','[64:2:76]'} ; % 10^-6 at 100ml/min
Input(48).BgConcentration = {'0','Log7','Log6'} ;

% vm7
% +/- background at several pulse concentrations
% small spikes and high rest potential
Input(49).cellname = '121212_1' ;
Input(49).OdorConcentration = '[10^-7,10^-6,10^-5,10^-4,10^-3]' ;
Input(49).OdorRsp{1} = {'[37:2:50,51:2:64]','[7:2:20,21:2:34]','[65:2:78]','[79:2:92]','[93:2:106]'} ; % no background
Input(49).OdorRsp{2} = {'[38:2:50]','[8:2:20]','[]','[]','[]'} ; % 10^-7 at 100ml/min
Input(49).OdorRsp{3} = {'[52:2:64]','[22:2:34]','[66:2:78]','[80:2:92]','[94:2:106]'} ; % 10^-6 at 100ml/min
Input(49).BgConcentration = {'0','Log7','Log6'} ;

% vm7
% +/- background at several pulse concentrations
% small spikes 
Input(50).cellname = '121213_1' ;
Input(50).OdorConcentration = '[10^-7,10^-6,10^-5,10^-4,10^-3]' ;
Input(50).OdorRsp{1} = {'[10:2:23]','[25:2:38]','[41:2:54]','[57:2:70]','[73:2:86]'} ; % no background
Input(50).OdorRsp{2} = {'[11:2:23]','[26:2:38]','[42:2:54]','[58:2:70]','[74:2:86]'} ; % 10^-6 at 100ml/min
Input(50).OdorRsp{3} = {'[]','[94:2:103]','[]','[]','[]'} ; % 10^-5 at 100ml/min (did not include control epochs (93:2:103) above)
Input(50).BgConcentration = {'0','Log6','Log5'} ;

% vm7
% +/- background at several pulse concentrations
% small spikes, several technical issues causing delay, several trials show firing rate oscillations
Input(51).cellname = '121220_1' ;
Input(51).OdorConcentration = '[10^-8,10^-7,10^-6,10^-5,10^-4,10^-3]' ;
Input(51).OdorRsp{1} = {'[10:2:23]','[26:2:39]','[41:2:54]','[69:2:79,81:2:88]','[91:2:104]','[107:2:120]'} ; % no background
Input(51).OdorRsp{2} = {'[11:2:23]','[27:2:39]','[42:2:54]','[70:2:79,82:2:88]','[92:2:104]','[108:2:120]'} ; % 10^-5 at 100ml/min
Input(51).OdorRsp{3} = {'[]','[]','[127:2:140]','[]','[]','[]'} ; % 10^-5 at 100ml/min (did not include control epochs (93:2:103) above)
Input(51).BgConcentration = {'0','Log5','wash'} ;

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
Input(53).OdorRsp{3} = {'[]','[]','[103:2:114]','[]','[]','[]'} ; % 10^-5 at 100ml/min (did not include control epochs (102:2:114) above)
Input(53).BgConcentration = {'0','Log5','wash'} ;

% vm7
% +/- background at several pulse concentrations
% good cell, short data blocks
Input(54).cellname = '121231_1' ;
Input(54).OdorConcentration = '[10^-8,10^-7,10^-6,10^-5,10^-4,10^-3]' ;
Input(54).OdorRsp{1} = {'[4:2:9]','[12:2:17]','[19:2:24]','[27:2:32]','[35:2:40]','[43:2:48]'} ; % no background
Input(54).OdorRsp{2} = {'[5:2:9]','[13:2:17]','[20:2:24]','[28:2:32]','[36:2:40]','[44:2:48]'} ; % 10^-5 at 100ml/min
Input(54).OdorRsp{3} = {'[]','[]','[52:2:56]','[]','[]','[]'} ; % 10^-5 at 100ml/min (did not include control epochs (51:2:56) above)
Input(54).BgConcentration = {'0','Log5','wash'} ;

% vm7
% +/- background at several pulse concentrations
% okay cell, small spikes, short data blocks
Input(55).cellname = '130108_1' ;
Input(55).OdorConcentration = '[10^-8,10^-7,10^-6,10^-5,10^-4,10^-3]' ;
Input(55).OdorRsp{1} = {'[4:2:9]','[12:2:17]','[19:2:24]','[27:2:32]','[35:2:40]','[43:2:48]'} ; % no background
Input(55).OdorRsp{2} = {'[5:2:9]','[13:2:17]','[20:2:24]','[28:2:32]','[36:2:40]','[44:2:48]'} ; % 10^-5 at 100ml/min
Input(55).OdorRsp{3} = {'[]','[]','[52:2:56]','[]','[]','[]'} ; % 10^-5 at 100ml/min (did not include control epochs (51:2:56) above)
Input(55).BgConcentration = {'0','Log5','wash'} ;



%% prep for export to hdf5 file and analyze

% prep matrix if it doesn't exist already
if ~exist('ForIgor')
    ForIgor=struct() ;
end

% analysis functions

% analysis of spontaneous activity in ORNS (data extraction from h5 file)
if perform.ORNanalysisA ;
    id = 'Spont' ;
    Temp = ORNanalysisA(Input,id,A)  ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.ORNanalysisB ;
    id = 'OdorRsp' ;
    Temp = ORNanalysisB(Input,id,A)  ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

% analysis of PN odor response in adapting long stim (I don't think this
% code exists yet!)
if perform.PNadaptationAnalysisA ;
    id = 'OdorRsp' ;
    Temp = PNadaptationAnalysisA ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end
    
% analysis of PN odor response +/- background
if perform.PNadaptationAnalysisB ;
    id1 = 'OdorRsp' ;
    id2 = 'OdorBgTimes' ;
    Temp = PNadaptationAnalysisB(Input,id1,id2,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.PNadaptationAnalysisC ;
    Temp = PNadaptationAnalysisC(Input,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.PNadaptationAnalysisD ;
    Temp = PNadaptationAnalysisD(Input,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end
 
if perform.PNadaptationAnalysisE ;
    Temp = PNadaptationAnalysisE(Input,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end
   
if perform.PNadaptationAnalysisF ;
    Temp = PNadaptationAnalysisF(Input,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.PNadaptationAnalysisG ;
    Temp = PNadaptationAnalysisG(Input,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.PNadaptationAnalysisH ;
    Temp = PNadaptationAnalysisH(Input,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

% analysis of odor response to random valve flicker
if perform.ValveNoiseAnalysis ;
    id = 'ValveNoise' ;
    Temp = ValveNoiseAnalysis(Input,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.ValveNoiseAnalysisB ;
    Temp = ValveNoiseAnalysisB(Input,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end
    
end % Data block (A) loop

%% across "cell" statistics
% if perform.PopDataFunctionXX ;
%     ForIgor = PopDataFunctionXX(ForIgor) ; % 
% end

if perform.PNadaptationAnalysisE ;
    Temp = PNadaptationAnalysisEpop(ForIgor) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.PNadaptationAnalysisH ;
    ForIgor = PopDataConcatAndStat(ForIgor,'average') ;
    ForIgor = PopDataConcatAndStat(ForIgor,'sem') ;
    Temp = PNadaptationAnalysisHpop(ForIgor) ;
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

if max(fieldLength)<32 ;
    cd Z:/cafaro' Documents'/Analysis/ForIgorHdfs/   
    exportStructToHDF5(ForIgor,'PNadaptationAnalysisH.h5','/')
else disp('name too long for igor')
end




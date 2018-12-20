% This script will save input info for cells being analyzed by different
% functions.

% Get the following outputs (1= get output, 0= don't) 
% assess cell spike response basic (CA only)
perform.CellCheckCA = 0 ;

% Burst and spike data from ca

perform.BStatCA = 0 ;

% interspike interval distribution
perform.ISIdisCA =0 ;

% variance of conductances
perform.VarGcheckExc = 0 ;
perform.VarGcheckInh = 0 ;

% Linear filter from light to conductance
perform.LinFExc = 0 ; % control G Exc
perform.LinFInh = 0 ; % control G Inh
perform.LinFInhApb = 0 ; % OFF channel inh isolated via APB

% Linear Nonlinear model from light to spike
perform.LNfiltersCA = 0 ; 

% Prespike stimuli from cell attached data
perform.STACA = 0 ;

% Prespike conductance (from Vclamp)
perform.STGCA = 0 ;

% prespike currents (from Iclamp)
perform.STIIC = 0 ;

% power spectrum of exc and inh
perform.PSg = 0 ;

% autocorrelation
perform.ACInh = 0 ;
perform.ACExc = 0 ;

% cross correlation
perform.CCexcinh = 0 ; % control exc and inh
perform.CCexcOffinh = 0 ; % control exc and apb inh 
perform.CCexcOninh = 0 ; % conrtrol exc and apb-control inh 

% alternating voltage experiments
perform.AltVanalysis = 0 ;
perform.AltVanalysisNegCorr = 0 ;
perform.AltVanalysisPosCorr = 0 ;

perform.AltVanalysisDS = 0 ;

perform.AltInhVanalysis = 0 ; % when inh is blocked
perform.AltInhVanalysisNegCorr = 0 ;
perform.AltInhVanalysisPosCorr = 0 ;

perform.AltVanalysis2 = 0 ;
perform.AltVanalysisNegCorr2 = 0 ;
perform.AltVanalysisPosCorr2 = 0 ;

perform.AltVanalysisDS3 = 1 ;

% IV steps
perform.IVSteps = 0 ;
perform.IVStepsAPV = 0 ;

% integrate model on simultaneous and shuffled conductances
perform.LIFforSimG = 0 ;

% dynamic clamp analysis on simultaneous conductances
perform.DCmidgetAnalyzer = 0 ;
perform.DCmidgetAnalyzer2 = 0 ;

perform.DCdsAnalyzer = 0 ;
perform.DCdsAnalyzer2 = 0 ;

% alternating voltage tau (time to get to exc rev pot)
perform.ValtTauAlt = 0 ;
perform.ValtTauAltDS = 0 ;
perform.ValtTauAltPosCorr = 0 ;
perform.ValtTauAltNegCorr = 0 ;

perform.ValtTauAltInh = 0 ;

% comparing control current step response with + drug
perform.icomparison = 0 ;

% pair cell attached data
perform.pairCADS = 0 ;

% pair whole cell data
perform.pairWCDSexcinh = 0 ;

% pair alternating voltage
perform.pairAltDS2 = 0 ;

% first order LN model light to G
perform.FirstOrderLNlight2G_Exc = 0 ;
perform.FirstOrderLNlight2G_Inh = 0   ;
perform.FirstOrderLNlight2G_InhApb = 0 ;

% dynamic clamp analysis of paired ds input +/-converging +/-pairwise correlations
perform.DCdsPairAnalyzer_141 = 0 ; 
perform.DCdsPairAnalyzer_157 = 0 ;
perform.DCdsPairAnalyzer_139 = 0 ;

perform.DCdsPairRepeatAnalyzer = 0 ;

% dynamic clamp analysis of paired ds input +/-converging +/-pairwise correlations
perform.DCdsPairAnalyzer3_141 = 0 ; 
perform.DCdsPairAnalyzer3_157 = 0 ; 
perform.DCdsPairAnalyzer3_139 = 0 ; 

% Parasol LIF model and G analysis for chapter 3 of thesis
perform.chapter3ParasolAnalysis = 0 ;
perform.chapter3ParasolDCAnalysis = 0 ;
perform.chapter3MidgetAnalysis = 0 ;

% using the following cells...
for A = [256:260] ;    
%% Cells for analysis  

Input(1).cellname = '060908Bc1' ; % took first set G for dc from this cell epochs(154:158 for exc) , Parasol with mulitple seeds repeated a few times each
Input(1).CA = '[79:108]' ;
Input(1).Exc = '[154:183]' ; % VC at inh reversal pot
Input(1).Inh = '[213:242]' ;
Input(1).InhApb = '[313:316]' ; % VC at exc revesal with APB present

Input(2).cellname = '060908Bc4' ; % 
Input(2).CA = '[45:74]' ;
Input(2).Exc = '[110:139]' ; % VC at inh reversal pot
Input(2).Inh = '[163:182]' ;
Input(2).InhApb = '[212:241]' ; % VC at exc revesal with APB present
Input(2).ExcApb = '[270:289]' ;
Input(2).InhWash = '[541:570]' ; % inh post apb
Input(2).ExcWash = '[616:645]' ;

Input(3).cellname = '061008Bc3' ; % 
Input(3).CA = '[87:126]' ;
Input(3).Exc = '[166:185]' ; % VC at inh reversal pot
Input(3).Inh = '[210:229]' ;
Input(3).InhApb = '[270:289]' ; % VC at exc revesal with APB present (30mV)
Input(3).InhApb2 = '[304:323]' ; % VC at another attempted reversal pot (40mV) 
Input(3).ExcApb = '[350:369]' ;

Input(4).cellname = '061008Bc4' ; % took second set of G for DC from this cell (epochs: 111:115,152:156,217:221), sfn poster example epochs(106,147,212) 
Input(4).CA = '[42:81]' ;
Input(4).Exc = '[106:125]' ; % VC at inh reversal pot
Input(4).Inh = '[147:166]' ;
Input(4).InhApb = '[212:251]' ; % VC at exc revesal with APB present
Input(4).ExcApb = '[292:325]' ;

Input(5).cellname = '061008Bc2' ; % 
Input(5).CA = '[25:64]' ;
Input(5).Exc = '[131:150]' ; % VC at inh reversal pot
Input(5).Inh = '[194:213]' ;

Input(6).cellname = '060908Bc2' ; % 
Input(6).CA = '[40:69]' ;
Input(6).Exc = '[109:138]' ; % VC at inh reversal pot
Input(6).Inh = '[165:174]' ;

Input(7).cellname = '060908Bc5' ; % 
Input(7).CA = '[39:68]' ;
Input(7).Exc = '[109:128]' ; % VC at inh reversal pot
Input(7).Inh = '[158:177]' ; % lost access and was not great throughout

Input(8).cellname = '060908Bc6' ; % 
Input(8).CA = '[30:59]' ; % did not have a blue mean on during stim
Input(8).Exc = '[91:120]' ; % VC at inh reversal pot
Input(8).Inh = '[184:193]' ;

Input(9).cellname = '061008Bc1' ; % 
Input(9).CA = '[50:89]' ;
Input(9).Exc = '[143:162]' ; % VC at inh reversal pot

Input(10).cellname = '060908Bc3' ; % 
Input(10).CA = '[45:74]' ;

% current clamped On parasol example
Input(11).cellname = '061908Bc1b' ;  % good patch 
Input(11).IC = '[34:43,110:124]' ;

Input(11).cellname = '030408Ac1' ; % midget cell from Fred
Input(11).Inh = '[116:135]' ;
Input(11).Exc = '[96:115]' ;

% Mouse cells alt Voltage experiments
Input(12).cellname = '020609Ec1' ; % ON-Off bistrat altvolt Rs = 18MOhms seems to be nice cell but possibly inaccurate reversal pots and short prestim time
Input(12).Steps = {'[222:226]','[217:221]','[212:216]','[227:231]','[232:236]','[237:241]','[242:246]','[247:251]'} ;
Input(12).SpepsPots = '[-80,-60,-45,-20,0,10,20,40]' ;
Input(12).Alt = '[180:194]' ; % attempted isolation (-50,15) reversal may have been closer to 5 or 10
Input(12).AltPosCorr = '[252:266]' ; % attempted positive correlation (beyond reversal pots) (-50,40)
Input(12).AltNegCorr = '[195:209]' ; % attempted negative correlation (between reversal pots) (-45,15)

Input(13).cellname = '020809Ec1' ; % on-off bistrat alt volt Rs=26 Mohms 
Input(13).Alt = '[69:83]' ; %(-50,5)
Input(13).AltNegCorr = '[134:148]' ; %(-40,5)
Input(13).Steps = {'[19:23]','[34:38]','[29:33,94:98]','[89:93]','[24:28,84:88]','[39:43]','[44:48]','[109:113]','[114:118]','[49:53,119:123]','[124:128]','[54:58,129:133]','[59:63]'} ;
Input(13).StepsPots = '[-60,-55,-50,-45,-40,-35,-25,-20,-10,-5,0,5,15]' ;

Input(14).cellname = '021209Ec1' ; % on-off bistrat alt volt Rs approx 20Mohms nice cell (alt v in apv) good example of gNMDA
Input(14).Steps = {'[85:89,165:169]','[160:164]','[80:84,155:159]','[180:184]','[150:154,175:179]','[170:174]','[90:94,145:149]','[140:144,185:194]','[95:99]','[135:139]','[130:134]','[100:104,125:129]','[115:124]','[110:114]','[105:109]'} ;
Input(14).StepsPots = '[-80,-70,-60,-55,-50,-45,-40,-30,-20,-10,-5,0,5,10,20]' ; % the holding potentials for the steps above
Input(14).StepsAPV = {'[260:264]','[265:269]','[270:274]','[275:279]','[365:369,375:379]','[280:284,355:359,370:374]','[360:364]','[250:259,285:289]','[290:294,350:354]','[345:349]','[295:299,330:334]','[300:304,335:339]','[305:309,340:344]','[310:314]','[315:319]','[320:324]','[325:329]'} ;
Input(14).StepsAPVPots = '[-80,-70,-60,-50,-45,-40,-30,-25,-20,-15,-10,-5,0,5,10,15,20]' ;
%Input(14).AltNegCorr = '[380:394]' ; %(-45,-20)
Input(14).AltNegCorr = '[425:439]' ; %(-45,-30)
Input(14).Alt = '[395:409,440:454]' ; %(-45,0)
Input(14).AltPosCorr = '[410:424]' ; %(-45,20)

Input(15).cellname = '021209Ec2' ; % on-off response but image may not be bistrat, seem to have lost access after steps causing contamination and correlations are negative 
Input(15).Steps = {'[20:24]','[15:19]','[30:34]','[25:29]','[35:39,55:59,80:84]','[40:44,50:54]','[45:49]','[60:64]','[65:69]','[70:74]','[75:79]'} ;
Input(15).StepsPots = '[-80,-60,-50,-40,-30,-20,-10,0,5,10,20]' ;
Input(15).StepsAPV = {'[135:139]','[140:144]','[145:149]','[150:154]','[160:164]','[155:159]','[130:134,165:169]','[170:174]','[175:179]','[180:184]','[185:189]','[190:194]','[195:199]','[200:204]'} ;
Input(15).StepsAPVPots = '[-80,-70,-60,-50,-45,-40,-30,-20,-10,-5,0,5,10,20]'
%Input(15).Alt = '[205:219]' ; % (-45,0) w/compensation as steps above
Input(15).Alt = '[230:244]' ; % (-45,0) no  compensation for all below too
%Input(15).Alt = '[328:342]' ; % (-45,0)
%Input(15).Alt = '[358:372]' ; % (-55,0) unclear if -45 or -55 is isolation
Input(15).AltPosCorr = '[245:259]' ; % (-45,20)
%Input(15).AltPosCorr = '[298:312]' ; % (-45,10)
%Input(15).AltPosCorr = '[373:387]' ; % (-55,20)
Input(15).AltNegCorr = '[260:274]' ; % (-45,-20) 
%Input(15).AltNegCorr = '[313:327]' ; % (-45,-10)
%Input(15).AltNegCorr = '[388:402]' ; % (-55,-20)

Input(16).cellname = '021609Ec1' ; % on-off bistrat, w/ strychnine (although not always?) bad access, several mistakes on this cell
Input(16).Steps = {'[115:119]','[110:114,120:124]','[125:129]','[130:134]','[135:139]','[140:144]','[145:149]','[150:154]','[155:159]','[160:164]','[165:169]'} ; % this previosly had APV ...
Input(16).StepsPots = '[-80,-60,-50,-40,-30,-20,-10,0,5,10,20]' ;
Input(16).StepsSty = {'[350:354]','[345:349,355:359]','[363:367]','[340:344]','[335:339]','[330:334]','[325:329]','[320:324]','[315:319,368:372]','[310:314,373:377,383:387]','[295:299,378:382]','[305:309]','[300:304]'} ; % these after 1min no perfusion
Input(16).StepsStyPots = '[-80,-60,-55,-50,-40,-30,-20,-10,0,10,20,30,40]' ;
Input(16).Alt = '[388:402]' ; %(-55,10)
Input(16).AltPosCorr = '[403:417]' ; %(-55,40)
Input(16).AltNegCorr = '[418:429]' ; %(-55,-20) lost access before recording 5 sets
Input(16).StepControl = '[165:169]' ; %(20)
Input(16).StepExp = '[295:299,378:382]' ; %(20 in strychnine)

Input(17).cellname = '020409Ec1' ; % on-off bistrat, 
% Input(17).Steps = {'[61:65]','[66:70]','[71:75]','[76:80]','[81:85]','[86:90]','[91:95]'} ; % 2.25Rh from dark
% Input(17).StepsPots = '[-60,-40,-20,0,10,20,30]' ;
Input(17).Steps = {'[106:115]','[116:125]','[126:135]','[136:145]','[146:155]'} ; % 2.25 from 2.25 mean
Input(17).StepsPots = '[-60,-40,-20,20,30]' ;
Input(17).Alt = '[160:171]' ; %(-60,30)
% Input(17)

Input(18).cellname = '020409Ec2' ; % on-off bistrat, 
Input(18).Steps = {'[64:68]','[74:78]','[69:73,134:138]','[59:63]','[54:58,139:143]','[144:148]','[94:98,149:153]','[84:88,174:178]','[79:83,154:158]','[89:93,169:173]','[159:163]'} ;
Input(18).StepsPots = '[-70,-65,-60,-50,-40,-20,0,10,20,30,40]' ;
Input(18).Alt = '[114:128]' ; %(-60,30)

Input(19).cellname = '020409Ec3' ; % on-off bistrat
Input(19).Steps = {'[53:57]','[48:52]','[58:62]','[83:87]','[78:82]','[73:77]'} ;
Input(19).StepsPots = '[-70,-60,-40,-20,0,20]' ;

Input(20).cellname = '030409Ec3' ; % on-off bistrat w/ gabazine
Input(20).AltInh = '[211:225]' ; % exc rev 40 
Input(20).AltPosCorrInh = '[266:280]' ; % hold 55 (true maybe 50)
Input(20).AltNegCorr = '[331:336]' ; % hold 35 (true maybe 50)

Input(21).cellname = '030409Ec2' ; % on-off bistrat
Input(21).Alt = '[116:130]' ;

Input(22).cellname = '030409Ec1' ; % off sustained
Input(22).Alt = '[128:142]' ; %(-40,20)
Input(22).AltPosCorr = '[239:253]' ; %(-60,20)

Input(23).cellname = '030609Ec2' ; % on-off bistrat + gabazine/stychnine and apv good cell
Input(23).AltInh = '[231:245]' ; % (-40,35)
Input(23).AltPosCorrInh = '[282:296]' ; %(-40,55)
Input(23).AltNegCorrInh = '[347:361]' ; % (-40,40)
Input(23).Steps = {'[70:74]','[75:79,90:94]','[85:89]','[80:84,95:99]','[100:104]','[105:114]','[115:124]'} ; % before gabazine,stychnine and apv,without a mean 
Input(23).StepsPots = '[-50,-40,-30,0,5,10,15]' ;
Input(23).Steps = {'[205:209]','[160:169]','[170:179]','[185:194]','[200:204]','[180:184,195:199]'} ; % after gabazine,stychnine and apv, without a mean
Input(23).StepsPots = '[-40,15,25,30,32,35]' ;
Input(23).StepControl = '[115:124]' ; %(15)
Input(23).StepExp = '[195:199]' ; %(15 in drugs)


% 030609Ec1 also has possibly useful data

Input(24).cellname = '031309Ec1' ; % on-off bistrat + gabazine/stychnine and apv good cell
Input(24).AltInh = '[402:416]' ; %(-50,25)
Input(24).AltPosCorrInh = '[453:467]'; %(-50,35)
Input(24).AltNegCorrInh = '[543:557]' ; %(-50,25) 
% Input(24).Steps = {'[60:69]','[70:79]','[80:84]','[85:89]','[90:94]','[95:99]'} ; % before gabazine/stychnine and apv, with a mean %THIS IS A DIFFERENT CELL!
% Input(24).StepsPots = '[-60,0,5,10,20,30]' ;
Input(24).Steps = {'[233:237]','[238,257]','[258:267]','[268:277]','[278:282]'} ; % before gabazine/stychinine and apv (although was present 15 minutes before), with a mean
Input(24).StepsPots = '[-60,0,10,20,30]' ;
Input(24).Steps = {'[363:377]','[348:352,378:387]','[353:362]','[338:347]',} ; % after inh blockade
Input(24).StepsPots = '[-50,20,25,30]' ;
Input(24).StepControl = '[268:277]' ; %(20)
Input(24).StepExp = '[348:352,378:387]' ; %(20 in drugs)

Input(25).cellname = '031309Ec2' ; % on-off bistrat + gabazine/stychnine and apv good cell
Input(25).AltInh = '[251:265]' ; %(-50,30)
%Input(25).AltInh = '[381:395]' ; %(-50,30)
Input(25).AltPosCorr = '[296:310]' ; %(-50,40)
Input(25).AltNegCorr = '[336:350]' ; %(-50,20)
Input(25).Steps = {'[65:74]','[85:89]','[90:99]','[75:84,101:104]','[105:109]','[110:114,120:124]','[115:119]'} ; % before gabazine/stychinine and apv, with a mean
Input(25).StepsPots = '[-60,-20,-10,0,10,20,30]' ;
Input(25).Steps = {'[210:214]','[170:179]','[240:244]','[180:194,245:249]','[195:209,225:239]'} ; % after inh blockade
Input(25).StepsPots = '[-50,20,25,30,35]' ;

Input(26).cellname = '031809Ec2' ; % on-off bistrat good cell (%some cell attached to noise stim)
Input(26).Steps = {'[72:76]','[57:61,87:91]','[202:206]','[82:86,182:186,197:201]','[192:196]','[62:66,77:81,187:191]','[67:71,167:171]','[162:166]','[152:156]','[92:96,147:151,157:161]','[97:101,142:146,172:176]','[102:106,177:181]','[107:111]','[112:121]'} ; 
Input(26).StepsPots = '[-80,-60,-55,-50,-45,-40,-20,-10,-5,0,10,20,30,40]' ;
Input(26).Alt = '[217:231]' ; %(-55,20)
Input(26).AltNegCorr = '[332:346]' ; %(-55,0)

Input(27).cellname = '041509Fc2b' ; % on-off bistrat in Ach block (this cell has step data too)
Input(27).Alt = '[30:44]' ;

Input(28).cellname = '041509Fc3' ; % on-off bistrat in Ach block (this cell has step data too)
Input(28).Alt = '[130:144]' ;

% below is midget alt v exp
Input(29).cellname = '042809Ec1' ; % very large on midget with big g, access may be limited by space clamp
Input(29).Alt = '[149:163]' ;

Input(30).cellname = '042809Ec4' ; % on midget, bad access
Input(30).Alt = '[184:198]' ;

Input(31).cellname = '042809Ec5' ; % on midget, bad access
Input(31).Alt = '[110:124]' ;

Input(32).cellname = '042909Ec1' ; % large on midget, good access but had on-off exc and inh
Input(32).Alt = '[72:86]' ;

Input(33).cellname = '042909Ec2' ; % off midget? unhealthy cell on-off g
Input(33).Alt = '[99:110]' ;

Input(34).cellname = '042909Ec3' ; % on midget, good access (rs =9)
Input(34).Alt = '[100:111]' ;

Input(35).cellname = '042909Ec4' ; % on midget + gabaznine and stychnine, bad access
Input(35).AltInh = '[211:225]' ;

Input(36).cellname = '042909Bc4' ; % on midget mulitclamp (Freds)
Input(36).ITC18flag = 1 ;
Input(36).Alt = '[62:76,151:165]' ;

Input(37).cellname = '042909Bc5' ; % On parasol multiclamp (rod light levels + correlations?)
Input(37).ITC18flag = 1 ;
Input(37).Alt = '[1066:1080]' ; 

Input(38).cellname = '060209Bc1' ; % on midget mulitclamp (very very weak light response - but good access)
Input(38).ITC18flag = 1 ;
Input(38).Alt = '[201:215,236:250]' ;

Input(39).cellname = '060809Ec2' ; % on midget multiclamp good access (good example)
Input(39).ITC18flag = 0 ;
%Input(39).Alt = '[196:225]' ; %(-50,15)
Input(39).Alt = '[246:281]' ; %(-50,20)
Input(39).Steps = {'[164:168]','[149:153]','[159:163]','[154:158]','[169:173]','[174:178]','[179:183]','[184:188]'} ;
Input(39).StepsPots = '[-70,-60,-50,-40,-10,0,10,15]' ;

Input(40).cellname = '060809Ec1' ; % on midget multiclamp very weak exc not great patch
Input(40).Alt = '[93:101]' ;
Input(40).Steps = {'[35:39]','[45:49]','[40:44]'} ;
Input(40).StepsPots = '[-60,-40,20]' ;

Input(41).cellname = '060809Ec3' ; % on midget (no image) multiclamp not great access
Input(41).Alt = '[117:131]' ; %(-55,20)
Input(41).Steps = {'[54:58]','[49:53]','[89:93]','[84:88]','[59:63]','[64:68]','[69:73]','[74:78]','[79:83]','[112:116]'} ;
Input(41).StepsPots = '[-80,-60,-55,-50,-40,-20,0,10,15,20]' ;

Input(42).cellname = '052609Bc1' ; % may be too big to be midget, access variable but not great
Input(42).ITC18flag = 1 ;
Input(42).Alt = '[87:101,148:162]' ;

Input(43).cellname = '061609Ec2' ; % midget (freds) multiclamp, 
Input(43).Alt = '[28:42]' ; %(-55,10) DIFFERNT SEEDS

Input(44).cellname = '061609Ec3' ; % on midget good access(freds) multiclamp, +/- gabazine and stychnine (this cell also has apv steps)
Input(44).Alt = '[31:45]' ; %(-50,10) DIFFERNT SEEDS
Input(44).AltNegCorr = '[46:60]' ; %(-40,10) DIFFERNT SEEDS
%Input(44).Alt = '[61:72]' ; %(-40,10) in INH block DIFFERNT SEEDS
%Input(44).Alt = '[74:85]' ; %(-40,30) in INH block DIFFERNT SEEDS

Input(45).cellname = '061609Ec4' ; % on midget good access (lost during inh block) (freds) +/- apv and +/- gabazine and strychnine, this cell has steps too, this cell has anticorrelated light driven mean g
%Input(45).Alt = '[20:34]' ; %(-55,5) control 
%Input(45).Alt = '[55:69]' ; %(-55,5) + APV
%Input(45).Alt = '[70:81]' ; %(-55,5) + INH block
Input(45).AltPosCorr = '[82:90]' ; %(-55,25) + INH block
%Input(45).AltPosCorr = '[91:105]' ; %(-55,45) + INH block

Input(46).cellname = '063009Ec2' ; % on midget multiclamp good access (13-17) + apv whole time
Input(46).Alt = '[116:130,186:200,216:230]' ; %(-50,20) APV
Input(46).AltNegCorr = '[146:160]' ; %(-30,20) APV
%Input(46).AltNegCorr = '[171:185]' ; %(-40,20) APV
Input(46).AltPosCorr = '[201:215]' ; %(-60,20) APV
Input(46).Steps = {'[16:20]','[21:25]','[26:30]','[31:35]','[36:40]','[41:45]'} ;
Input(46).StepsPots = '[-80,-60,-40,-20,0,20]' ;
Input(46).StepsAPV = {'[76:80]','[81:85]','[86:90]','[71:75]','[91:95]','[106:110]','[101:105]','[96:100]'} ;
Input(46).StepsAPVPots = '[-80,-60,-40,-20,0,10,15,20]' ;

Input(47).cellname = '063009Ec3' ; % on midget multiclamp good access (13-18) but not as good during apv
Input(47).Alt = '[86:115]' ; %(-50,20) control
%Input(47).Alt = '[157:171]' ; %(-50,20) apv

Input(48).cellname = '063009Ec4' ; % on midget multiclamp good access (12-15) +/- apv and +/- sychnine and gabazine
Input(48).Alt = '[80:94]' ; %(-50,20) control
Input(48).AltNegCorr = '[99:113]' ; %(-40,20)
%Input(48).Alt = '[159:173,189:203,219:227]' ; %(-50,20) +APV
%Input(48).AltNegCorr = '[174:188]' ; %(-40,20) +APV
Input(48).AltPosCorr = '[204:218]' ; %(-60,20) +APV
Input(48).AltInh = '[276:290]' ; %(-50,20) + INH block
Input(48).StepControl = '[236:240]' ; %(20,+apv)
Input(48).StepExp = '[266:270]' ; %(20 +apv, stychnine and gabazine)

Input(49).cellname = '070709Ec1' ; % on midget multiclamp good access +/- strychnine and gabazine
Input(49).Alt = '[134:145]' ; %(-60,20)
%Input(49).Alt = '[147:161]' ; %(-50,20)
Input(49).AltInh = '[264:278]' ; %(-50,20) +INH block
%Input(49).AltPosCorrInh = '[295:309]' ; %(-50,35) +INH block
Input(49).AltPosCorrInh = '[317:331]' ; %(-50,45) +INH block 
Input(49).StepControl = '[167:171]' ; %(20)
Input(49).StepExp = '[222:226]' ; %(20 stychnine and gabazine)

Input(50).cellname = '070709Ec2' ; % on midget multiclamp good access +/- strychnine and gabazine
Input(50).Alt = '[70:78,81:86]' ; %(-50,20)
Input(50).AltInh = '[149:157]' ; %(-50,20) +INH block
Input(50).AltPosCorrInh = '[159:173]' ; %(-50,40) + INH block
Input(50).StepControl = '[94:98]' ; %(20)
Input(50).StepExp = '[124:128]' ; %(20 stychnine and gabazine)

Input(51).cellname = '072809Ec1' ; % midget (off?) dynamic clamp access not good and no drug block, WAVEFORMS CUT SHORT!
Input(51).ITC18flag = 0 ; % indicates interval
Input(51).dcSimShuffle = '[5:20]' ; 

Input(52).cellname = '080609Bc1' ; % on-off ds cell low mean light but good access and light response (STIMULUS TIMING VERY INPRECISE!)
Input(52).ITC18flag = 1 ;
Input(52).AltDS = '[88:102]' ; % no mean bar around 4000 p*/mcone (-50,10)
Input(52).AltDS = '[103:132]' ; % mean around 800 P*/mcone (-50,10)
Input(52).AltDS = '[133:162]' ; % mean around 800 (-50,20)

Input(53).cellname = '081209Bc1' ; % on-off ds cell, access and light response okay, low mean light
Input(53).ITC18flag = 1 ;
Input(53).CADS = '[43:47]' ; % cell attached
Input(53).AltDS = '[78:107]' ;%(-50,15)
%Input(53).AltDS = '[109:123]' ; %(-50,15)
%Input(53).AltDS = '[125:154]' ; %(-50,20)


Input(54).cellname = '063009Bc4b' ; % off midget alt-v (freds)
Input(54).ITC18flag = 1 ;
Input(54).Alt = '[0:14]' ; %(-45,20)

Input(55).cellname = '063009Bc5' ; % off midget alt-v (freds)
Input(55).ITC18flag = 1 ;
Input(55).Alt = '[12:17]' ; %(-55,10) mean 1V
Input(55).Alt = '[18:26]' ; %(-55,10) mean .4V

Input(56).cellname = '072809Bc2' ; % midget dynamic clamp (freds), decent access, synaptic block, naked retina, WAVEFORMS CUT SHORT!
Input(56).ITC18flag = 1 ; % indicates interval
Input(56).dcSimShuffle = '[21:36,37:52]' ;

Input(57).cellname = '081809Ec1' ; % on midget (sparse), decent access
Input(57).dcSimShuffle = '[94:125]' ;

Input(58).cellname = '081809Ec2' ; % on midget (sparse), decent access
Input(58).dcSimShuffle = '[32:79]' ; 

Input(59).cellname = '081809Ec3' ; %on midget, decent access, used rs compensation
Input(59).dcSimShuffle = '[31:46]' ;

Input(60).cellname = '081809Ec4' ; %on midget, decent access
Input(60).dcSimShuffle = '[50:65]' ;

Input(61).cellname = '082509Bc1' ; %on midget, not good access
Input(61).ITC18flag = 1 ;
Input(61).AltNegCorr = '[78:86]' ; %(-50,20)
Input(61).Alt = '[118:132]' ; %(-50,30)
Input(61).AltNegCorr = '[135:149,155:169]' ; %(-35,30)

Input(62).cellname = '082509Bc2' ; %on midget probably but sparse, poor access
Input(62).ITC18flag = 1 ;
Input(62).Alt = '[55:69]' ; %(-50,20)
Input(62).AltNegCorr = '[73:87]' ; %(-30,20)

Input(63).cellname = '082609Bc1' ; %on midget not great access
Input(63).ITC18flag = 1 ;
Input(63).CAnoise = '[10:19] ' ; %cell attached noise (other cells may have similar data)
Input(63).Alt = '[40:54]' ; %(-50,20)
Input(63).AltNegCorr = '[70:84]' ; %(-35:20)
Input(63).AltNegCorr = '[85:99]' ; %(-40:20)

Input(64).cellname = '082609Bc2' ; %on midget not great access
Input(64).ITC18flag = 1 ;
Input(64).CAnoise = '[10:19]' ; %cell attached noise
Input(64).Alt = '[45:59]' ; %(-60,20)
%Input(64).Alt = '[70:84]' ; %(-60,30)
Input(64).AltNegCorr = '[105:119]' ; %(-30,30)
%Input(64).AltNegCorr = '[120:134]' ; %(-40,30)

Input(65).cellname = '082609Bc3' ; %on midget good access +/- gabazine and strychnine
Input(65).ITC18flag = 1 ;
%Input(65).Alt = '[39:53]' ; %(-60,15)
Input(65).Alt = '[73:87] ' ; %(-60,20)
Input(65).AltInh = '[118:132]' ; %(-60,20) +INH BLOCK
Input(65).AltPossCorrInh = '[159:173] ' ; %(-60,55) +INH BLOCK
Input(65).StepControl = '[88:92]' ; %(20)
Input(65).StepExp = '[108:112]' ; %(20 stychnine and gabazine)

Input(66).cellname = '082609Ec1' ; %on midget (freds) good access - this cell has longer stim time than others
%Input(66).Alt = '[26:40]' ; %(-50,5)
Input(66).Alt = '[78:92]' ; %(-60,5)
Input(66).AltPosCorr = '[115:123]' ; %(-60,15)

Input(67).cellname = '082609Ec2' ; %on midget (freds) good access +/- gabazine and stychnine
Input(67).Alt = '[13:27]' ; %(-50,10)
Input(67).AltNegCorr = '[30:44]' ; %(-40,10)
Input(67).AltInh = '[80:94]' ; %(-50,10) INH BLOCK
Input(67).AltPosCorrInh = '[101:115]' ; %(-50,16) INH BLOCK
Input(67).AltPosCorrInh = '[116:130]' ; %(-50,26) INH BLOCK
Input(67).ControlStep = '[45:49]' ; %inh step before inh block at rev used
Input(67).DrugStep = '[70:74]' ; % inh step after inh block
Input(67).StepControl = '[45:49]' ; %(10)
Input(67).StepExp = '[70:74]' ; %(10 stychnine and gabazine)

Input(68).cellname = '090809Bc1' ; %on midget (probably but weird spike response and intrinsic props), good access
Input(68).ITC18flag = 1 ;
Input(68).dcSimShuffle = '[44:75]' ; 

Input(69).cellname = '090809Bc2' ; %on midget (probably but weird spike response and intrinsic props), good access
Input(69).ITC18flag = 1 ;
Input(69).dcSimShuffle = '[46:77]' ; 

Input(70).cellname = '091509Bc2' ; % on midget (big and not def on midget), okay access, this cell has inward current in presence of NBQX+APV
Input(70).ITC18flag = 1 ;
Input(70).Alt = '[50:64]' ; %(-60,25)
Input(70).AltNegCorr = '[75:89]' ; %(-40,25)

Input(71).cellname = '100709Bc1' ; % on midget alt-v dynamic clamp with adjusted exp design (use DCmidgetAnalyzer2-for all beyond), okay access some change in firing rate
Input(71).ITC18flag = 1 ;
Input(71).dcSimShuffle = '[58:81]' ;
Input(71).dcSimRepeat = '[58,82:85]' ; 

Input(72).cellname = '101309Bc2' ; % on midget alt-v dc, good access
Input(72).ITC18flag = 1 ;
Input(72).dcSimShuffle = '[30:53]' ;
Input(72).dcSimRepeat = '[30,54:59]' ;

Input(73).cellname = '101309Ec3' ; % on midget (Freds) alt-v dc, "decent"
Input(73).dcSimShuffle = '[2:25]' ; % 10,20
Input(73).dcSimRepeat = '[2,26:34]' ; 
% Input(73).dcSimShuffle = '[35:58]' ; 
% Input(73).dcSimRepeat = '[35,59:66]' ;
% Input(73).dcSimShuffle = '[67:90]' ;
% Input(73).dcSimRepeat = '[67,91:100]' ;
% Input(73).dcSimShuffle = '[101:124]' ; % 15,30
% Input(73).dcSimRepeat = '[101,125:131]' ;
% Input(73).dcSimShuffle = '[132:155]' ;
% Input(73).dcSimRepeat = '[132,156:165]' ;

Input(74).cellname = '101309Ec2' ; % on midget (Freds) "lousy"
Input(74).dcSimShuffle = '[4:27]' ;
Input(74).dcSimRepeat = '[4,28]' ;

Input(75).cellname = '101409Ec2' ; % on midget (Freds) "nice"
Input(75).dcSimShuffle = '[1:24]' ; % 15,30
Input(75).dcSimRepeat = '[1,25:33]' ; 
% Input(75).dcSimShuffle = '[34:57]' ; 
% Input(75).dcSimRepeat = '[34,58:64]' ;
% Input(75).dcSimShuffle = '[65:88]' ; 
% Input(75).dcSimRepeat = '[65,89:94]' ;
% Input(75).dcSimShuffle = '[95:118]' ;
% Input(75).dcSimRepeat = '[95,119:125]' ;

Input(76).cellname = '101409Bc1' ; % on midget access variable
Input(76).ITC18flag = 1 ;
Input(76).dcSimShuffle = '[0:23]' ; % 10,20
Input(76).dcSimRepeat = '[0,24:26]' ; 
% Input(76).dcSimShuffle = '[38:61]' ; % less stable than first set
% Input(76).dcSimRepeat = '[38,62:65]' ;

Input(77).cellname = '101309Bc1' ; % on midget access not great, weak dc response
Input(77).ITC18flag = 1 ;
Input(77).CAnoise = '[5:14]' ; 
Input(77).dcSimShuffle = '[35:58]' ; %no drugs
Input(77).dcSimRepeat = '[35,59:63]' ;
% Input(77).dcSimShuffle = '[83:106]' ; %blockers
% Input(77).dcSimRepeat = '[83,107:111]' ;

Input(78).cellname = '092309Bc3' ; % on-off ds cell decent access
Input(78).ITC18flag = 1 ;
Input(78).CADS = '[0:4]' ;
Input(78).AltDS = '[47:76]' ; %(-55,30)
%Input(78).AltDS = '[-45,30]' ;%(-45,30)

Input(79).cellname = '102009Bc1' ; % on midget alt-v dc, good access but weak spike rates
Input(79).ITC18flag = 1 ;
Input(79).dcSimShuffle = '[53:76]' ; %(10,20)
Input(79).dcSimShuffle = '[79:102]' ; %(15,30)
Input(79).dcSimRepeat = '[79,103:108]' ;
% Input(79).dcSimShuffle = '[113:136]' ; %(15,30)
% Input(79).dcSimRepeat = '[113,137:146]' ;
% Input(79).dcSimShuffle = '[147:170]' ; %(15,30)
% Input(79).dcSimRepeat = '[171:174]' ;

Input(80).cellname = '102009Ec1' ; % on midget (Freds) good access
Input(80).dcSimShuffle = '[2:25]' ; %(15,30)
Input(80).dcSimRepeat = '[2,26:31]' ; 
% Input(80).dcSimShuffle = '[32:55]' ;
% Input(80).dcSimRepeat = '[32,56:64]' ;
% Input(80).dcSimShuffle = '[65:88]' ;
% Input(80).dcSimRepeat = '[65, 89:98] ' ;

Input(81).cellname = '102009Ec2' ; % on midget (Freds) "ok"
Input(81).dcSimShuffle = '[14:37] ' ; %(12,24)
Input(81).dcSimRepeat = '[14,38:45]' ; 

Input(82).cellname = '102109Fc1' ; % on midget good access good example
Input(82).CAnoise = '[8:17]' ;
Input(82).dcSimShuffle = '[34:57]' ; %(10,20)
Input(82).dcSimRepeat = '[34,58:63]' ;
% Input(82).dcSimShuffle = '[64:87]' ;
% Input(82).dcSimRepeat = '[64,88:93]' ;
% Input(82).dcSimShuffle = '[94:117]' ;
% Input(82).dcSimRepeat = '[118:123]' ;

Input(83).cellname = '102609Bc1' ; % on parasol access not great (this cell actually on the 27th)
Input(83).ITC18flag = 1 ;
Input(83).SfCA = '[15:1145]' % cell attached split field stim, cell drifted and was moved
Input(83).SfExc = '[1212:1256,1457:1606]' ;
Input(83).SfInh = '[1177:1211,1307:1456]' ;

Input(84).cellname = '102709Bc2' ; % on parasol access not great
Input(84).ITC18flag = 1 ;
Input(84).SfCA = '[21:866]' ; 
Input(84).SfExc = '[867:876,942:1006,1160:1309]' ;
Input(84).SfInh = '[892:941,1010:1159]' ;

Input(85).cellname = '110409Bc1' ; % mouse on-off ds cell good access moving bar (used for dynamic clamp exp)
Input(85).ITC18flag = 1 ;
%Input(85).CADS = '[23:27]' ; % mean .1 V cell 
Input(85).CADS = '[28:37]' ; % mean .2 V
Input(85).AltDS = '[79:99,100:114,115:138]' ;%(-60,20) varied Rs , mean .2V
%Input(85).AltDS = '[139:171]' ;% zero mean

Input(86).cellname = '110409Bc2' ; % on-off ds cell good access
Input(86).ITC18flag = 1 ;
Input(86).CADS = '[40:49]' ; % mean .2 V
%Input(86).CADS = '[50:54]' ;% mean 0
Input(86).AltDS = '[85:117]' %,152:184,185:196]' ; % mean .2V
%Input(86).AltDS = '[118:150]' ;

Input(87).cellname = '110409Bc3' ; % on-off ds cell good access
Input(87).ITC18flag = 1 ;
Input(87).CADS = '[2:11]' ; % .2 mean PULSE AT SAME TIME (MISTAKE)
%Input(87).CADS = '[12:16]' ; % 0 mean PULSE AT SAME TIME (MISTAKE)
Input(87).AltDS = '[47:91,92:109]' ; %.2 V mean

Input(88).cellname = '111109Bc3' ; % on-off ds moving bar good access
Input(88).ITC18flag = 1 ;
Input(88).dcSimShuffle = '[65:104]' ; %(10 20) no drugs
Input(88).dcSimRepeat = '[65,65]' ; % no repeat data attained on this cell
%Input(88).dcSimShuffle = '[129:138]' ; %(10 20) no drugs

Input(89).cellname = '111109Bc4' ; % on-off ds moving bar good access
Input(89).ITC18flag = 1 ;
Input(89).dcSimShuffle = '[30:69]' ; %(7,14) no drugs
Input(89).dcSimRepeat = '[30,70:74]' ;

Input(90).cellname = '111109Bc3' ; % on-off ds moving bar bad access
Input(90).ITC18flag = 1 ;
Input(90).dcSimShuffle = '[40:49]' ;
Input(90).dcSimRepeat = '[40,40]' ; % no repeat data attained on this cell

Input(91).cellname = '111109Bc1' ; % on-off ds moving bar bad access
Input(91).ITC18flag = 1 ;
Input(91).dcSimShuffle = '[14:23]' ;
Input(91).dcSimRepeat = '[14,14]' ; % no repeat data attained on this cell

Input(92).cellname = '112509Bc1' ; % on-off ds moving bar good access (looked for spikelets in this cell)
Input(92).ITC18flag = 1 ;
Input(92).IclampDsDc = '[207:211]' ; % w/ mean light current clamp moving bar
Input(92).CADS = '[10:159]' ; % epochs to assess spike count distribution in dark (Watch out for break ins!) 
Input(92).dcSimShuffle = '[160:199]' ;
Input(92).dcSimRepeat = '[160,200:205]' ; %(5 15)
%Input(92).dcSimShuffle = '[226:245]' ; %(5,20)

Input(93).cellname = '112509Bc2' ; % on-off ds moving bar good access (looked for spikelets in this cell)
Input(93).ITC18flag = 1 ;
Input(93).CADS = '[10:14,25:29]' ; % no mean
Input(93).CADS = '[15:24]' ; % mean .2V (2/3 bar)
Input(93).dcSimShuffle = '[46:85]' ; %(3,7)
Input(93).dcSimRepeat = '[46,86:94]' ; %
%Input(93).dcSimShuffle = '[95:114]' ; %(Ginh reversal = -60mV)
%Input(93).dcSimShuffle = '[115:134]' ; %(held closer to spike threshold)

Input(94).cellname = '112509Bc3' ; % on-off ds moving bar good access (good example) 
Input(94).ITC18flag = 1 ;
Input(94).CADS = '[38:42]' ; % no mean
Input(94).CADS = '[43:47]' ; % mean .2V (2/3 bar)
%Input(94).dcSimShuffle = '[50:69]' ; %(4,8)
Input(94).dcSimShuffle = '[78:117]' ; %(4,12 used steady current inj instead of steady hold slow variable inhj as above)
Input(94).dcSimRepeat = '[118:125]' ; % 4 12

Input(95).cellname = '112509Ec1' ; % on-off (Freds) "good access"
Input(95).dcSimShuffle = '[30:69]' ; %(5,15)
Input(95).dcSimRepeat = '[30,70:76]' ;
%Input(95).dcSimShuffle = '[77:106]' ;

Input(96).cellname = '112509Ec2' ; % on-off (Freds) "decent" 
Input(96).dcSimShuffle = '[5:44]' ; %(7,21)
Input(96).dcSimRepeat = '[5,45:49]' ; 

Input(97).cellname = '112509Ec3' ; % on-off? (Freds) - had spikelet and unclear on-off bistrat!!!!
Input(97).dcSimShuffle = '[15:34]' ; %(5,15)
Input(97).dcSimRepeat = '[15,35:38]' ;

Input(98).cellname = '112509Ec4' ; % on-off (Freds)
Input(98).dcSimShuffle = '[48:87]' ; %(3,9)
Input(98).dcSimRepeat = '[48,88:92,133:136]' ; 
Input(98).dcSimShuffle = '[93:132]' ; %(3,9)

% on-off ds pairs below
Input(99).cellname = '022610Bc2' ; % on-off ds pair moving bar, good responses, about 50uM apart no image
Input(99).ITC18flag = 1 ;
Input(99).step = '[97:101]' ; % step w/ background
Input(99).CADS = '[102:114,116:119,121:145,147:148]' ; % moving bar pair cell attached

Input(100).cellname = '030210Bc1' ; % on-off ds pair moving bar, good responses, neighbors, see image
Input(100).ITC18flag = 1 ;
Input(100).step = '[47:51]' ; % step w/o background
Input(100).step = '[52:56]' ; % step w background
Input(100).CADS = '[57:64,66:95]' ; 

Input(101).cellname = '031610Fc1' ; % (Freds) on-off ds pair moving bar,no background mean  
Input(101).CADS = '[2:16]' ;
Input(101).ExcDS = '[22:27]' ; 

Input(102).cellname = '031610Fc1b' ; % same cell as 101
Input(102).InhDS = '[3:12]' ;

Input(103).cellname = '031610Fc2' ; % (Freds) on-off ds pair moving bar
Input(103).CADS = '[11:45,50:66]' ;

Input(104).cellname = '031610Fc3' ; % (freds) on-off ds pair moving bar (one cell not DS!)
Input(104).CADS = '[8:27]' ;
Input(104).ExcDS = '[28:32]' ;

Input(105).cellname = '040610Bc1' ; % pair on-off ds good access, see image
Input(105).ITC18flag = 1 ;
Input(105).CADS = '[10:24]' ;
%Input(105).InhDS = '[33:42]' ; % 10,10
Input(105).InhDS = '[66:75]' ; % 20,20
Input(105).ExcDS = '[43:52,86:96]' ; % -60
%Input(105).ExcDS = '[97:98]' ; % -40
%Input(105).InhExcDS = '[53:55]' ; % 10,-60
Input(105).InhExcDS = '[56:65,100:109]' ; % 20,-60
Input(105).ExcInhDS = '[76:85,110:119]' ; % -60,20

Input(106).cellname = '041610Bc1' ; % pair on-off ds ca, see image
Input(106).ITC18flag = 1 ;
Input(106).CADS = '[129:143,144:203]' ;
%Input(106).CADS = '[144:203]' ;

Input(121).cellname = '032310Fc1' ; % (Freds) pair on-off ds, no image
Input(121).CADS = '[5:25]' ;

Input(107).cellname = '032310Fc2' ; % (Freds) pair on-off ds, bad access one cell, no image
Input(107).CADS = '[4:24]' ;
Input(107).ExcDS = '[26:34]' ; 

Input(108).cellname = '040610Bc2' ; % ds pair, no image
Input(108).CADS = '[7:26]' ;

Input(109).cellname = '041610Fc1' ; % (Freds) pair ds, one cell a little wierd, no image
Input(109).CADS = '[6:14,16:25]' ;
Input(109).AltVdsPair = '[31:56,72:89]' ; %(-50,10/15)
%Input(109).AltVdsPair = '[57:71]' ; %(15:50)

Input(110).cellname = '041610Fc2' ; % ds pair (Freds)
% ???

Input(111).cellname = '041610Fc2b' ; % cont from 041610Fc2
% ???

Input(112).cellname = '041610Fc2c' ; % cont from 041610Fc2, no image
Input(112).AltVdsPair = '[2:31]' ;
Input(112).AltVdsPair = '[32:45]' ;

Input(113).cellname = '042010Bc2' ; % ds pair, good access both but ca rundown, see image
Input(113).ITC18flag = 1 ;
Input(113).CADS = '[20:52]' ;
Input(113).AltVdsPair = '[73:108]' ; %(-55,10)

Input(114).cellname = '042110Bc1' ; % ds pair, ca watch out for 2nd cell same amp, wc only one cell with access, see image
Input(114).ITC18flag = 1 ;
Input(114).CADS = '[50:89]' ; % no mean?
Input(114).AltDS = '[103:126]' ;

Input(115).cellname = '042110Bc2' ; % ds pair, ca watch for rundown, good access on 2nd cell but may be pair useable. see image
Input(115).ITC18flag = 1 ;
Input(115).CADS = '[10:26,28:49]' ;
Input(115).AltVdsPair = '[65:94]' ;


% midget cell attached data from Fred for psth for sup figure 5 of my first paper
Input(116).cellname = '030408Ac1' ;
Input(116).amp = 0 ;
Input(116).CA = '[5:33]' ; % amp 1

Input(117).cellname = '030608Ac1' ;
Input(117).amp = 1 ;
Input(117).CA = '[27:41]' ; % amp 2

Input(118).cellname = '030608Ac2' ;
Input(118).amp = 0 ;
Input(118).CA = '[19:28]' ;  % amp 1

Input(119).cellname = '032008Ac1' ; 
Input(119).amp = 0 ;
Input(119).CA = '[12:40]' ; % amp 1
% end of Freds group of ca midget data

Input(120).cellname = '042010Fc1c' ; % (Freds) on-off ds pair, some ca also but lost notes (042010Fc1), no image
Input(120).AltVdsPair = '[9:32]' ; %(-50,25)
Input(120).AltVdsPair = '[44:61]' ; %(35,-50)
%Input(120).AltVdsPair = '[62:88]' ; %(-50,35)

% 121 is above
Input(122).cellname = '082710Fc1' ; % on-off ds pair (90 apart), ca only, many from a few directions
Input(122).CADS = '[63:383,386:412]' ;

Input(123).cellname = '082710Fc2' ; % on-off ds pair (180 apart), ca only, see crappy image
Input(123).CADS = '[8:207]' ; % epochs with different background light

Input(124).cellname = '090910Fc1' ; % on-off ds pair, ca and alt-v, good data, see image
Input(124).CADS = '[15:170]' ; % 
%Input(124).AltVdsPair = '[215:370]' ; %(-50,10)
Input(124).AltVdsPair = '[371:430]' ; %(-50,20)
%Input(124).AltVdsPair = '[233:307,311:340,344:370]' ; %(-50,10) 

Input(125).cellname = '090910Fc2' ; % on-off ds pair, ca , good data, see image
Input(125).CADS = '[13:170]' ; % 

% parasol data with repeated seed
Input(126).cellname = '010908Bc1' ; % parasol in APV
Input(126).Exc = '[65:69,85:89]' ; 
Input(126).Inh = '[70:74,90:99]' ;

Input(127).cellname = '100610Bc1' ; % pair on-off ds ca and alt-v good data, cell 1 not great access, see image
Input(127).ITC18flag = 1 ;
Input(127).CADS = '[17:158]' ; 
Input(127).AltVdsPair = '[173:271]' ; %(-60,20)
Input(127).AltVdsPair = '[272:328]' ; %(-60,20)

Input(128).cellname = '102610Bc1' ; % on-off ds ca and alt-v at around 1rh* nackground, single cell but needs to be analyzed as pair, good access but maybe not compensated
Input(128).ITC18flag = 1 ;
Input(128).CADS = '[19:28]' ;
Input(128).AltVdsPair = '[50:97]' ;

Input(129).cellname = '102610Bc2' ; % on-off ds ca and alt-v at around 1rh* background, single cell but needs to be analyzed as pair, good access
Input(129).ITC18flag = 1 ;
Input(129).CADS = '[10:18]' ; % this has pair data
Input(129).AltVdsPair = '[35:121]' ; 

Input(130).cellname = '110310Fc3' ; % on-off ds pair ca alt-v at around 1rh* background, Freds
Input(130).CADS = '[32:38]' ;
Input(130).AltVdsPair = '[42:92]' ;
Input(130).ExcDS = '[93:103]' ; 
Input(131).InhDS = '[104:119]' ;

Input(131).cellname = '110310Fc2' ; % on-off ds pair ca alt-v at around 1rh*, Freds, not much data

Input(132).cellname = '102910Bc1' ; % on-off ds pair ca, alt-v around 1rh*, bad access on both and weak light responses
Input(132).ITC18flag = 1 ;
Input(132).CADS = '[6:15]' ; 
Input(132).AltVdsPair = '[44:73]' ;

Input(133).cellname = '110210Bc1' ; % on-off ds ca and alt-v at around 1rh* background, single cell but needs to be analyzed as pair (cell 2 only), good access
Input(133).ITCflag = 1 ;
Input(133).CADS = '[7:16]' ; 
Input(133).AltVdsPair = '[62:109]' ;

Input(134).cellname = '102910Bc1' ; % on-off ds ca and alt-v around 1rh* background around, bad access on both and weak tunging curves, somas about 30um apart
Input(134).ITC18flag = 1 ;
Input(134).CADS = '[6:15]' ; % ERROR LED PULSE MAY HAVE BEEN ON 
Input(134).AltVdsPair = '[44:73]' ; 

Input(135).cellname = '110510Bc1' ; % on-off ds ca and alt-v around 1rh* background 10x bar, single cell but analyze in pair, good access but bad light response,

Input(136).cellname = '110910Bc1' ;  % on-off ds ca and alt-v around 1rh* background, single cell wc but analyze in pair, good access,
Input(136).ITC18flag = 1 ;
Input(136).CADS = '[13:22]' ; 
Input(136).AltVdsPair = '[62:106]' ; % cell 1 only

Input(137).cellname = '111010Bc1' ; % on-off ds ca and alt-v around 50Rh* background 100 contrast, single cell but analyze as pair, good access
Input(137).ITC18flag = 1 ;
Input(137).CADS = '[13:22]' ;
Input(137).AltDS = '[43:69]' ; % watch for lost access
Input(137).channel = 0 ;

Input(138).cellname = '111210Fc1' ; % on-off ds pair (Freds) ca, ***
Input(138).CADS = '[24:33]' ;

Input(139).cellname = '111210Fc1b' ; % same pair as above ,alt-v
Input(139).AltVdsPair = '[10:90]' ;

Input(140).cellname = '111210Fc1c' ; % same pair as above , v-clamp gj assessment
Input(140).VstepPair = '[0:74]' ;

Input(141).cellname = '111710Fc1' ; % on-off ds pair (Freds), alt-v
Input(141).CADS = '[12:16,17:26]' ; % both cells may not have been simultaneously recorded
Input(141).AltVdsPair = '[37:39,43:54,58:63,67:90,94:96]' ; % look for issues!!!!!! bars are in sets of 2

Input(142).cellname = '111710Fc1b' ; % same pair as above, v-clamp gj assessment (def coupled)
Input(142).VstepPair = '[0:209]' ;

Input(143).cellname = '111910Fc1' ; % on-off ds pair (Freds), alt-v, good image
Input(143).CADS = '[16:20]' ;
Input(143).AltVdsPair = '[25:27,31:54]' ; % notes that something is wrong
Input(143).InhDS = '[68:101]' ;

Input(144).cellname = '111910Fc1b' ; % same pair as above, gj assessment with cnqx and strychnine
Input(144).VstepPair = '[0:89]' ;
Input(144).VstepPair = '[94:218]' ; % +cnqx and stychnine

Input(145).cellname = '120310Fc2' ; % on-off ds pair, ca and alt-v good data (Freds)
Input(145).CADS = '[0:3]' ;
Input(145).AltVdsPair = '[10:27,31:39]' ; % messed up single hold trial included

Input(146).cellname = '120310Fc2b' ; % same pair as above
Input(146).AltVdsPair = '[0:17,21:47]' ;

Input(147).cellname = '120310Fc2c' ; % same pair as above 
Input(147).AltVdsPair = '???' ; % has some 8 direction 

Input(148).cellname = '120310Fc2d' ; % same pair as above
Input(148).VstepPair = '[10:99]' ;

Input(150).cellname = '120310Bc1' ; % on-off ds single alt v, (recorded in pair)
Input(150).ITC18flag = 1 ;
Input(150).AltDS = '[49:123]' ; 
Input(150).channel = 0 ;

Input(151).cellname = '120810Fc1' ; % on-off ds single alt v, (recorded as pair - early data on pair), not good access on either
Input(151).CADS = '[10:14]' ; % cell 1 only on amp 1
Input(151).CADS = '[20:24]' ; % cell 2 only on amp 1
Input(151).AltVdsPair = '[42:53]' ; % pair
Input(151).AltDS = '[55:78]' ; % okay access on cell 2 only
Input(151).channel = 1 ;

Input(152).cellname = '120810Fc2' ; % on-off ds pair alt-v good access but synaptic input dropped off a cliff at some point
Input(152).CADS = '[13:17]' ; % cell 1 amp 1
Input(152).CADS = '[18:26]' ; % cell 2 amp 2
Input(152).AltVdsPair = '[48:86]' ; % more epochs but exc in cell 1 only stable until epoch 86 

Input(153).cellname = '120810Fc2b' ; % same pair as above
Input(153).VstepPair = '[1:81]' ;

Input(154).cellname = '121010Fc2' ; % on-off ds alt-v good access (Fred and I), computer issue spread data in several files
Input(154).CADS = '[17:21]' ; 
Input(154).AltVdsPair = '[27:32]' ; % 6 repeats

Input(155).cellname = '121010Fc2b' ; % same pair as above
Input(155).AltVdsPair = '[1:21]' ; % 21 repeats

Input(156).cellname = '121010Fc2e' ; % same pair as above
Input(156).AltVdsPair = '[1:9]' ; % 9 repeats

Input(157).cellname = '121610Bc1' ; % on-off ds alt-v decent access (tissue not very good though),
Input(157).ITC18flag = 1 ;
Input(157).CADS = '[10:14]' ; % cell2 amp 2
Input(157).CADS = '[20:24]' ; % cell1 amp 2
Input(157).AltVdsPair = '[46:90,94:108]' ; %'[46:90,94:120]' ; % seemed to drop response strength in cell 1 after epoch 108

Input(158).cellname = '121610Bc1b' ; % same pair as above
Input(158).ITC18flag = 1 ;
Input(158).VstepPair = '[12:71]' ;

Input(159).cellname = '020411Bc1' ; % on-off; pair ds dc, good access, long recording
Input(159).ITC18flag = 1 ;
%Input(159).dc141 = '[26:44]' ; % incomplete (8/8,16/16)
Input(159).dc141 = '[49:172]' ; % complete (8/8,16/16)
%Input(159).dc141 = '[205:286]' ; % incomplete (8/8,16/16)
%Input(159).dc141 = '[288:395]' ; % incomplete + synaptic blockers (8/8,16/16)
Input(159).dcRepeat = '[49:130,205:286]' ; %141

Input(160).cellname = '012711Fc1' ; % on-off; pair ds dc, poor access, look for complete loss in last couple epochs
Input(160).dc141 = '[23:98]' ; % incomplete

Input(161).cellname = '011111Bc1' ; % on-off ds; pair ds dc, good access, note different exc inh ratios, looked for spikelets - none seen, DATA ACTUALLY COLLECTED 2/11/11 not 1/11/11 
Input(161).ITC18flag = 1 ;
Input(161).dc157 = '[38:225]' ; % (8/8,16/16) complete
%Input(161).dc157 = '[353:540]' ; % (8/8,24/24) complete
Input(161).dc141 = '[229:352]' ; % (8/8,24/24) complete

Input(162).cellname = '011111Bc2' ; % on-off ds; pair ds dc, good access, looked for spikelets - none seen, DATA ACTUALLY COLLECTED 2/11/11 not 1/11/11 
Input(162).ITC18flag = 1 ;
Input(162).dc139 = '[17:120]' ; % complete (8/8,16/16)
%Input(162).dc139 = '[121:224]' ; % complete (8/8,16/16)
Input(162).dcRepeat = '[17:224]' ; % 139 

Input(163).cellname = '011111Bc3' ; % on-off ds; pair ds dc, good access, note different inh rev pots and exc inh ratios, looked for spikelets - none seen, DATA ACTUALLY COLLECTED 2/11/11 not 1/11/11 
Input(163).ITC18flag = 1 ;
Input(163).dc141 = '[21:144]' ; % complete (8/8,16/16)
%Input(163).dc141 = '[147:270]' ; % complete (8/8,16/16) % inh rev at -80 (all previous ds pair dc data using -60mV
%Input(163).dc141 = '[271:394]' ; % complete (12/12,16/16) % inh rev at -80

Input(164).cellname = '020911Bc1' ; % on-off (ds?); pair ds dc, access okay, spiklets during dc?
Input(164).ITC18flag = 1 ;
Input(164).dc157 = '[13:44]' ; % incomplete (8/8,16/16)
%Input(164).dc157 = '[45:116]' ; % incomplete (8/8,16/16)
Input(164).dcRepeat = '[13:116]' ; % 157

Input(165).cellname = '020911Bc2' ; % unknown (prob not on-off ds);  pair ds dc, good access, spiklets seen during hyperpol. light steps
Input(165).ITC18flag = 1 ;
Input(165).dc157 = '[17:204]' ; % complete (8/8,16/16)
%Input(165).dc157 = '[205:352]' ; % complete (8/8, 16/16)
Input(165).dcRepeat = '[17:352]' ; % 157

Input(166).cellname = '021711Bc1' ; % on-off (?); pair ds dc, okay-good access, FROM HERE OUT SHUFFLE IS DIFFERENT IN IGOR CODE
Input(166).ITC18flag = 1 ;
Input(166).dc141 = '[18:145]' ; % complete (8/8,16/16)
%Input(166).dc141 = '[146:273]' ; % complete (8/8,16/16)
Input(166).dc139 = '[284:387]' ; % complete (8/8,16/16)
%Input(166).dc139 = '[388:491]' ; % complete (8/8,16/16)
Input(166).dcRepeat = '[18:273]' ; % 141
Input(166).dcRepeat = '[284:491]' ; % 139

Input(167).cellname = '021711Bc2' ; % on-off ds;pair ds dc, okay access, 
Input(167).ITC18flag = 1 ;
Input(167).dc141 = '[275:402]' ; % complete (15/15,16/16)
%Input(167).dc141 = '[403:530]' ; % complete (15/15,16/16)
%Input(167).dc141 = '[531:658]' ; % complete (15/15,16/16)+synaptic blockers
Input(167).dcRepeat = '[275:530]' ; % 141

Input(168).cellname = '021811Bc1' ; % poor spiking

Input(169).cellname = '031111Bc1' ; % on-off ds; pair ds dc, good access but this cell showed evidence for spikelets.
Input(169).ITC18flag = 1 ;
%Input(169).dc141 = '[12:139]' ; % complete (8/8,16/16)
Input(169).dc141 = '[140:267]' ; % complete (12/12,16/16)
%Input(169).dc141 = '[268:395]' ; % complete (12/12,16/16)
Input(169).dcRepeat = '[140:395]' ; % 141

Input(170).cellname = '031111Bc2' ; % on-off ds; pair ds dc, good access but this cell may have showed evidence for spikelets.
Input(170).ITC18flag = 1 ;
Input(170).dc141 = '[21:148]' ; % complete (12/12,16/16)
%Input(170).dc141 = '[149:276]' ; % complete (12/12,16/16)
%Input(170).dc141 = '[282:409]' ; % complete (12/12,16/16)
Input(170).dcRepeat = '[21:276]' ; % 141
Input(170).dcRepeat = '[149:276,282:409]' ; % 141

Input(171).cellname = '031111Bc3' ; % on-off ds; pair ds dc, okay access
Input(171).ITC18flag = 1 ;
%Input(171).dc141 = '[108:235]' ; % complete (12/12,16/16)
Input(171).dc141 = '[238:365]' ; % complete (16/16,16/16)
%Input(171).dc141 = '[366:493]' ; % complete but tissue ran out of solution at some unknown epoch (16/16,16/16)
Input(171).dcRepeat = '[238:493]' ; % 141

Input(172).cellname = '031011Bc1' ; % on-off ds; pair ds dc, poor access
Input(172).ITC18flag = 1 ;
Input(172).dc141 = '[53:124]' ; % incomplete

Input(173).cellname = '031711Bc1' ; % on-off ds; pair ds dc, good access but cell and spike threshold began to hyperpolarize, so stoped recording
Input(173).ITC18flag = 1 ;
Input(173).dc139 = '[25:96]' ; % (8/8,16/16)incomplete

Input(174).cellname = '031711Bc1' ; % on-off ds; pair ds dc, okay access, looked for spikelets at the end
Input(174).ITC18flag = 1 ;
Input(174).dc139 = '[131:234]' ; %(12/12,16/16) complete
%Input(174).dc139 = '[235:338]' ; %(12/12,16/16) complete
%Input(174).dc139 = '[339:442]' ; %(12/12,16/16) complete
Input(174).dcRepeat = '[131:338]' ; % 139
Input(174).dcRepeat = '[235:442]' ; % 139

Input(175).cellname = '031711Bc2' ; % on-off ds; pair ds dc, good access but lots of spikelets during dc
Input(175).ITC18flag = 1 ;
Input(175).dc139 = '[89:192]' ; %(12/12,16/16) complete?-BUT bad gain setting on last 2 epochs. AMP CRASHED AT THE END SO NOT SURE IF LAST COUPLE EPOCHS WERE RECORDED

Input(176).cellname = '033011Bc1' ; % on-off ds; pair ds dc, okay access, no spikelets seen
Input(176).ITC18flag = 1 ;
Input(176).dc139 = '[15:118]' ; % (8/8,16/16) complete
%Input(176).dc139 = '[121:224]' ; % (8/8,16/16) complete

Input(177).cellname = '033011Bc2' ; % on-off ds; pair ds dc, okay access but poor spike response
Input(177).ITC18flag = 1 ;
Input(177).dc139 = '[23:126]' ; % (10/10,20/20) complete
%Input(177).dc139 = '[127:230]' ; % (10/10,20/20) complete
%Input(177).dc139 = '[231:334]' ; % (10/10,15/15) complete

Input(178).cellname = '033011Bc3' ; % on-off ds; pair ds dc, good access but poor spike response, had to use very high e-i ratio
Input(178).ITC18flag = 1 ;
%Input(178).dc139 = '[61:164]' ; % (10/10,8/8) complete
Input(178).dc139 = '[195:298]' ; % (20/20,5/5) complete
%Input(178).dc139 = '[299:402]' ; % (20/20,5/5) complete

Input(179).cellname = '033111Bc1' ; % on-off ds?; pair ds dc, good access but poor spiking 
Input(179).ITC18flag = 1 ;
Input(179).dc139 = '[18:121]' ; % (10/10,10/10) complete

Input(180).cellname = '033111Bc2' ; % on-off ds; pair ds dc, okay access, good spiking 
Input(180).ITC18flag = 1 ;
Input(180).dc139 = '[28:131]' ; % (10/10,15/15) complete
%Input(180).dc139 = '[132:235]' ; % (10/10,15/15) complete

Input(181).cellname = '033111Bc3' ; % on-off ds; pair ds dc, good access, okay spiking , I was not present at the end
Input(181).ITC18flag = 1 ;
%Input(181).dc157 = '[20:100]' ; % (10/10,12/12) complete
Input(181).dc157 = '[105:296]' ; % (10/10,10/10) complete
%Input(181).dc157 = '[297:377]' ; % (10/10,10/10) complete
Input(181).dc139 = '[381:484]' ; % (15/15,15/15) complete

Input(182).cellname = '040611Bc1' ; % on-off ds; pair ds dc, okay access, okay spiking, spiking fell off in second set
Input(182).ITC18flag = 1 ;
Input(182).dc157 = '[23:214]' ; % (10/10,15/15) complete
%Input(182).dc157 = '[215:326]' ; % (10/10,15/15) incomplete

Input(183).cellname = '040611Bc2' ; % on-off ds; pair ds dc, okay access, good spiking, spiking fell off in second set 
Input(183).ITC18flag = 1 ;
Input(183).dc157 = '[15:206]' ; % (10/10,15/15) complete
%Input(183).dc157 = '[212:243]' ; % (10/10,15/15) incomplete

Input(184).cellname = '040611Bc3' ; % on-off ds; pair ds dc, good access, okay spiking, looked for spikelets briefly at the end 
Input(184).ITC18flag = 1 ; 
Input(184).dc157 = '[15:206]' ; % (10/10, 15/15) complete
%Input(184).dc157 = '[207:262]' ; % (10/10, 15/15) incomplete

Input(185).cellname = '040611Bc4' ; % on-off ds; exc/ inh +/- gabazine during light steps

Input(186).cellname = '041411Bc1' ; % on-off ds pair, alt-v during moving bar, decent access
Input(186).ITC18flag = 1 ;
Input(186).AltVdsPair = '[43:102]' ;

Input(187).cellname = '050611Fc1' ; % on-off ds pair, alt-v, cell 2 lost access between epochs 64-96, cell 1 could not be compensated, see image
Input(187).ITC18flag = 1 ;
Input(187).AltVdsPair = '[49:96]' ;

Input(188).cellname = '050611Fc2' ; % on-off ds, alt-v, recorded as pair but only cell 2 was alive, okay access
Input(188).ITC18flag = 1 ;
Input(188).AltDS = '[32:94]' ; % look at channel 2 only
Input(188).channel = 1 ;

Input(189).cellname = '051111Fc1b' ; % on-off ds pair, alt-v, poor access cell 1, okay in cell 2, lost notes for ca tuning curves, see image
Input(189).ITC18flag = 1 ;
Input(189).AltDS = '[10:93]' ; 
Input(189).channel = 1 ;

Input(190).cellname = '051711Fc1' ; % on-off ds pair (Freds), alt-v, "access okay" nice light responses
Input(190).ITC18flag = 1 ;
Input(190).AltDS = '[26:94,97:126]' ; 
Input(190).channel = 0 ; 

Input(191).cellname = '051711Bc1b' ; % same pair as above, not coupled
Input(191).ITC18flag = 1 ;
Input(191).VstepPair = '[0:119]' ;

% post fixing conductances for on-off g clamp experiments
Input(192).cellname = '072211Bc1' ; % on-off ds dc, good access
Input(192).ITC18flag = 1 ;
Input(192).gHWflag = 1 ;
Input(192).dc141 = '[29:156]' ; %(8,8/16,16) complete
Input(192).dc139 = '[161:264]' ; %(8,8/16,16) complete
%Input(192).dc157 = '[265:348]' ; %(8,8/16,16) incomplete poor spiking
%Input(192).dc157 = '[350:425]' ; %(10,10/12,12) incomplete poor spiking
%Input(192).dc157 = '[427:442]' ; %(25,25/20,20) incomplete poor spiking

Input(193).cellname = '080311Bc1' ; % on-off ds dc, yfp+ cell, good access, good spiking
Input(193).ITC18flag = 1 ;
Input(193).gHWflag = 1 ;
Input(193).dc157 = '[7:182]' ; %(10,10/18,18) complete 

Input(194).cellname = '080311Bc2' ; % on-off ds dc, yfp+ cell, good access, good spiking
Input(194).ITC18flag = 1 ;
Input(194).gHWflag = 1 ;
Input(194).dc141 = '[17:144]' ; %(10,10/20,20) complete
Input(194).dc139 = '[160:263]' ; %(10,10/15,15) complete

Input(195).cellname = '080311Bc3' ; % on-off ds dc, yfp+ cell, good access, good spiking
Input(195).ITC18flag = 1 ;
Input(195).gHWflag = 1 ;
%Input(195).dc141 = '[211:338]' ; %(10,10/20,20) complete but POOR SPIKING AND ACCESS
Input(195).dc157 = '[33:208]' ; %(10,10/15,15) complete 

Input(196).cellname = '080311Bc4' ;  % on-off ds dc, yfp+ cell, good access, good spiking
Input(196).ITC18flag = 1 ;
Input(196).gHWflag = 1 ;
%Input(196).dc141 = '[296:390]' ; %(10,10/20,20) incomplete, CELL CRASHED DURING THIS SET but not sure when
Input(196).dc139 = '[191:294]' ; %(10,10/20,20) complete
Input(196).dc157 = '[13:188]' ; %(10,10/20,20) complete 

Input(197).cellname = '080311Bc5' ; % on-off ds dc, yfp+ cell, okay access, okay spiking
Input(197).ITC18flag = 1 ;
Input(197).gHWflag = 1 ;
Input(197).dc141 = '[26:153]' ; %(10,10/15,15) complete

Input(198).cellname = '080311Bc6' ; % on-off ds dc, yfp+ cell, good access, good spiking
Input(198).ITC18flag = 1 ;
Input(198).gHWflag = 1 ;
Input(198).dc141 = '[287:377]' ; %(10,10/20,20) incomplete and spiking died on this block
Input(198).dc139 = '[6:109]' ; %(10,10/20,20) complete
Input(198).dc157 = '[111:286]' ; %(10,10/20,20) complete 

Input(199).cellname = '080311Bc7' ; % on-off ds dc, yfp+ cell, good access, good spiking
Input(199).ITC18flag = 1 ;
Input(199).gHWflag = 1 ;
Input(199).dc141 = '[180:307]' ; %(10,10/20,20) complete
Input(199).dc157 = '[4:179]' ; %(10,10/20,20) complete 

% PAIRS RECORDED on 080811 have significant timing jitter from bad gain on
% frame monitor.  MUST align by a consistent upsweep.  
Input(200).cellname = '080811Bc1' ; % on-off ds pair altV moving bar, different dirrections but not sure if 90 or 180 (only 1 is yfp+)
Input(200).ITC18flag = 1 ;
Input(200).triggerMiss = ones(1,75)*2 ; % number of missed upsweeps (can be vector if different for each epoch)
%Input(200).triggerMiss([21,42,96,102]) = 1 ;
Input(200).triggerMiss([21,42]) = 1 ;
%Input(200).AltVdsPair = '[29:133]' ; % (-60,20) lost access cell 1 before cell 2 but could be good until 103
Input(200).AltVdsPair = '[29:103]';

Input(201).cellname = '080811Bc2' ; % on-off ds pair altV moving bar, same direction both yfp+. good access
Input(201).ITC18flag = 1 ;
Input(201).triggerMiss = 1 ;
Input(201).AltVdsPair = '[90:158]' ; % (-50,20) lost access cell 2 before cell 2 after 104 at some time

Input(202).cellname = '080811Bc3' ; % same pair as above, they don't look connected
Input(202).ITC18flag = 1 ;
Input(202).triggerMiss = 1 ;
Input(202).VstepPair = '[0:25,27:126]' ;

Input(203).cellname = '080811Bc4' ; % on-off ds pair altV moving bar, same direction both yfp+. not great access, no background
Input(203).ITC18flag = 1 ;
Input(203).triggerMiss = 2 ;
%Input(203).AltVdsPair = '[64:78,80:181]' ; % (-50,20) lost access cell 2 before cell 2 after 104 at some time
Input(203).AltVdsPair = '[110:169]' ; % these epochs look like better access in cell 1

Input(204).cellname = '080811Bc5' ; % same pair as above, they don't look connected
Input(204).ITC18flag = 1 ;
Input(204).triggerMiss = 1 ;
Input(204).VstepPair = '[1:110]' ;

% single cell from same day below

Input(205).cellname = '080811Bc6' ; % on-off ds altV, single cell, recorded as pair, moving bar, good access
Input(205).ITC18flag = 1 ;
Input(205).triggerMiss = 2 ;
Input(205).AltDS = '[41:115]' ; % (-50,20)
Input(205).channel = 0 ;

Input(206).cellname = '081111Bc2' ; % on-off ds dc, old stuff repeated with fixed g
Input(206).ITC18flag = 1 ;
Input(206).gHWflag = 1 ;
Input(206).dcSimShuffle = '[45:84]' ; % (10,10)
Input(206).dcSimRepeat = '[45:47]' ; % not really repeated data

Input(207).cellname = '080811Bc1' ; % on-off ds pair altV moving bar, different dirrections but not sure if 90 or 180 (only 1 is yfp+)
Input(207).ITC18flag = 1 ;
Input(207).triggerMiss = ones(1,75)*2 ; % number of missed upsweeps (can be vector if different for each epoch)
Input(207).triggerMiss([21,42]) = 1 ;
Input(207).AltDS = '[29:103]' ; % (-60,20) lost access cell 1 before cell 2 but could be good until 103
Input(207).channel = 0 ;

Input(208).cellname = '080811Bc1' ; % on-off ds pair altV moving bar, different dirrections but not sure if 90 or 180 (only 1 is yfp+)
Input(208).ITC18flag = 1 ;
Input(208).triggerMiss = ones(1,75)*2 ; % number of missed upsweeps (can be vector if different for each epoch)
Input(208).triggerMiss([21,42]) = 1 ;
Input(208).AltDS = '[29:103]' ; % (-60,20) lost access cell 1 before cell 2 but could be good until 103
Input(208).channel = 1 ;

Input(209).cellname = '080811Bc2' ; % on-off ds pair altV moving bar, same direction both yfp+. good access
Input(209).ITC18flag = 1 ;
Input(209).triggerMiss = 1 ;
Input(209).AltDS = '[90:158]' ; % (-50,20) lost access cell 2 before cell 2 after 104 at some time
Input(209).channel = 0 ;

Input(210).cellname = '080811Bc2' ; % on-off ds pair altV moving bar, same direction both yfp+. good access
Input(210).ITC18flag = 1 ;
Input(210).triggerMiss = 1 ;
Input(210).AltDS = '[90:158]' ; % (-50,20) lost access cell 2 before cell 2 after 104 at some time
Input(210).channel = 1 ;

Input(211).cellname = '080811Bc4' ; % on-off ds pair altV moving bar, same direction both yfp+. not great access, no background
Input(211).ITC18flag = 1 ;
Input(211).triggerMiss = 2 ;
Input(211).AltDS = '[110:169]' ; % (-50,20) lost access cell 2 before cell 2 after 104 at some time
Input(211).channel = 0 ;

Input(212).cellname = '080811Bc4' ; % on-off ds pair altV moving bar, same direction both yfp+. not great access, no background
Input(212).ITC18flag = 1 ;
Input(212).triggerMiss = 2 ;
Input(212).AltDS = '[110:169]' ; % (-50,20) lost access cell 2 before cell 2 after 104 at some time
Input(212).channel = 1 ;

Input(213).cellname = '081911Bc1b' ; % yfp+ likely on-off ds, dc old stuff redo, okay access, spiking often outside of injected g
Input(213).ITC18flag = 1 ;
Input(213).gHWflag = 1 ;
Input(213).dcSimShuffle = '[13:52]' ; % (10,10)
Input(213).dcSimRepeat = '[13:15]' ; % not really repeated data

Input(214).cellname = '081911Bc5' ; % yfp+ likely on-off ds, dc old stuff redo, decent access but drift in spiking and then cell died
Input(214).ITC18flag = 1 ; 
Input(214).gHWflag = 1 ;
Input(214).dcSimShuffle = '[87:118]' ; % incomplete up to trial 16
Input(214).dcSimRepeat = '[87:89]' ; % not really repeated data

Input(215).cellname = '082411Bc1b' ; % on-off ds, dc old-stuff redo, okay access, some little drift in rest, cell died quickly, incomplete
Input(215).ITC18flag = 1 ; 
Input(215).gHWflag = 1 ;
Input(215).dcSimShuffle = '[8:33]' ; % (14,20) incomplete up to trial 13
Input(215).dcSimRepeat = '[8:9]' ; % not really repeated data

Input(216).cellname = '082411Bc2' ; % on-off ds, dc old-stuff redo, good access but weak spiking
Input(216).ITC18flag = 1 ; 
Input(216).gHWflag = 1 ;
Input(216).dcSimShuffle = '[21:60]' ; % (15,20) 
Input(216).dcSimShuffle = '[66:105]' ; % (18,20) better tuning but worse access 
Input(216).dcSimRepeat = '[21:22]' ; % not really repeated data

Input(217).cellname = '082411Bc3' ; % on-off ds, dc old-stuff redo, good access but weak spiking 
Input(217).ITC18flag = 1 ; 
Input(217).gHWflag = 1 ;
Input(217).dcSimShuffle = '[13:52]' ; % (15,20) 
Input(217).dcSimShuffle = '[56:95]' ; % (25,20) maybe better tuning 
Input(217).dcSimRepeat = '[13:14]' ; % not really repeated data

Input(218).cellname = '082411Bc3b' ; % I clamp looking for spikelets during hyperpolarization

Input(219).cellname = '082411Bc4' ; % on-off ds, dc old-stuff redo, good access but lost it at the end
Input(219).ITC18flag = 1 ; 
Input(219).gHWflag = 1 ;
Input(219).dcSimShuffle = '[14:53]' ; % (13,20) 
Input(219).dcSimRepeat = '[14:15]' ; % not really repeated data

Input(220).cellname = '091511Bc3' ; % on-off ds pair altV moving bar, Freds - "good access, nice responses"
Input(220).ITC18flag = 1 ;
Input(220).triggerMiss = 1 ;
Input(220).CADS = '[26:30]' ; % cell 1
Input(220).CADS = '[31:35]' ; % cell 2
%Input(220).AltVdsPair = '[42:71]'; %(-55,20) stimulus issues?
Input(220).AltVdsPair = '[84:113]' ; %(-55,30)
% Input(220).AltVdsPair = '[114:128]' ; %(-55,40)

Input(221).cellname = '091511Bc2' ; % on-off ds pair altV moving bar, Freds "right next to each other, access ok but bot great"
Input(221).ITC18flag = 1 ;
Input(221).triggerMiss = 1 ;
Input(221).CADS = '[0:4]' ; % cell 1
Input(221).CADS = '[5:9]' ; % cell 2
Input(221).AltVdsPair = '[16:75]' ; %(-55,25)

% older on parasol dynamic clamp data (222-235 data blocks below)
Input(222).cellname = '070908Bc1' ; % ON parasol good access 
Input(222).dcControl = '[261:265,305:309,310:314,315:319]' ; %(much more data out of blockers or different sample rate)
Input(222).dcMinusInh = '[266:270,275:279]' ;

Input(223).cellname = '092408Ac2' ; %ON parasol okay access - not useful data

Input(224).cellname = '092408Ac3' ; %ON parasol okay access - not usedul data

Input(225).cellname = '092408Bc1' ; %ON parasol good access
Input(225).dcControl = '[32:36,52:56,57:61,67:71,82:86,102:106,149:168]' ;
Input(225).dcMinusInh = '[37:51,62:66,87:91,107:111]' ;
Input(225).dcMinusExc = '[]' ;

Input(226).cellname = '092408Bc2' ; %ON parasol variable access
Input(226).dcControl = '[50:54,60:64,70:74,95:99,100:104]' ;
Input(226).dcMinusInh = '[55:59,65:69]' ;

Input(227).cellname = '092408Bc3' ; %ON parasol good access
Input(227).dcControl = '[27:31,32:36,47:51]' ;
Input(227).dcMinusInh = '[17:21,22:26,37:41,42:46,52:56]' ;

Input(228).cellname = '102108Bc1' ; % ON parasol not great access DIFFERENT SET OF G USED BELOW

Input(229).cellname = '102108Bc2' ; %ON parasol okay access
Input(229).dcControl = '[78:82,108:112,128:132,153:157]' ;
Input(229).dcMinusInh = '[88:92,113:117,133:137,158:162]' ;

Input(230).cellname = '102108Bc3' ; %ON parasol okay access
Input(230).dcControl = '[63:67,68:72,88:92,93:97,113:117,138:142,143:147]' ;
Input(230).dcMinusInh = '[73:77,98:102,123:127]' ;

Input(231).cellname = '102108Bc4' ; %ON parasol okay access

Input(232).cellname = '102108Bc5' ; %ON parasol okay access

Input(233).cellname = '102108Bc6' ; %ON parasol okay access

Input(234).cellname = '102208Bc1' ; %ON parasol okay access

Input(235).cellname = '102208Bc2' ; %ON parasol okay access

Input(236).cellname = '092211Bc1' ; % (Freds) decent ds pair "screwed up light level"
Input(236).ITC18flag = 1 ;
Input(236).triggerMiss = 0 ;
Input(236).CADS = '[8:12]' ; % cell 1
Input(236).CADS = '[4:7]' ; % cell 2
Input(236).AltVdsPair = '[30:74,77:112]' ; %(-55,20)

Input(237).cellname = '092211Bc2b' ; % (Freds) good ds pair, (vstep data also)
Input(237).ITC18flag = 1 ;
Input(237).triggerMiss = 1 ;
Input(237).CADS = '[0:3]' ; % cell 1
Input(237).CADS = '[4]' ; % cell 2 ????
Input(237).AltVdsPair = '[8:142]' ; %(-55,20)

Input(238).cellname = '092011Bc1' ; % Freds "decent ON parasol" Gratings LOOK FOR LOST ACCESS
Input(238).ITC18flag = 1 ;
Input(238).CA = '[0:37]' ;
Input(238).Exc = '[70:79,100:109,115:144]' ;
Input(238).Inh = '[80:99,145:154]' ;

Input(239).cellname = '092011Bc2' ; % Freds "nice ON parasol" Gratings 
Input(239).ITC18flag = 1 ;
Input(239).CA = '[10:62]' ;
Input(239).Exc = '[96:141,251:300,332:366]' ;
Input(239).Inh = '[81:95,154:198,301:331]' ;

Input(240).cellname = '092011Bc3' ; % Freds "Off parasol - ok" Gratings 
Input(240).ITC18flag = 1 ;
Input(240).CA = '[0:29]' ;
Input(240).Exc = '[]' ;
Input(240).Inh = '[]' ;

Input(241).cellname = '092011Bc4' ; % Freds "nice Off parasol" Gratings +/- ly/apb 
Input(241).ITC18flag = 1 ;
Input(241).CA = '[]' ;
Input(241).Exc = '[]' ;
Input(241).Inh = '[]' ;

Input(242).cellname = '092011Bc5' ; % Freds "nice On parasol" Gratings +/- ly/apb 

% single cell ds
Input(243).cellname = '092211Bc2b' ; % (Freds) good ds pair
Input(243).ITC18flag = 1 ;
Input(243).triggerMiss = 1 ;
Input(243).channel = 0 ; % cell 1
Input(243).AltDS = '[8:142]' ; %(-55,20)

Input(244).cellname = '092211Bc2b' ; % (Freds) good ds pair
Input(244).ITC18flag = 1 ;
Input(244).triggerMiss = 1 ;
Input(244).channel = 1 ; % cell 1
Input(244).AltDS = '[8:142]' ; %(-55,20)

Input(245).cellname = '091511Bc3' ; % on-off ds pair altV moving bar, Freds - "good access, nice responses"
Input(245).ITC18flag = 1 ;
Input(245).triggerMiss = 1 ;
Input(245).channel = 0 ;
Input(245).AltDS = '[84:113]' ; %(-55,30)

Input(246).cellname = '091511Bc3' ; % on-off ds pair altV moving bar, Freds - "good access, nice responses"
Input(246).ITC18flag = 1 ;
Input(246).triggerMiss = 1 ;
Input(246).channel = 1 ;
Input(246).AltDS = '[84:113]' ; %(-55,30)

Input(247).cellname = '091511Bc2' ; % on-off ds pair altV moving bar, Freds "right next to each other, access ok but bot great"
Input(247).ITC18flag = 1 ;
Input(247).triggerMiss = 1 ;
Input(247).channel = 0 ;
Input(247).AltDS = '[16:75]' ; %(-55,25)

Input(248).cellname = '091511Bc2' ; % on-off ds pair altV moving bar, Freds "right next to each other, access ok but bot great"
Input(248).ITC18flag = 1 ;
Input(248).triggerMiss = 1 ;
Input(248).channel = 1 ;
Input(248).AltDS = '[16:75]' ; %(-55,25)

% pair data
Input(249).cellname = '100711Bc1' ; % on-off ds pair alt-v moving bar, Fred, 1 +yfp 1-yfp, "nice ds pair gain of second cell 20 for constant hold at end- may be messed up in file (first cell gain 10)
Input(249).ITC18flag = 1 ;
Input(249).triggerMiss = 1 ;
Input(249).CADS = '[0:3]' ; % cell 1
Input(249).CADS = '[4:7]' ; % cell 2?
Input(249).AltVdsPair = '[12:86]' ; 

Input(250).cellname = '100711Bc1b' ; % pair above v-step data
Input(250).ITC18flag = 1 ;
Input(250).VstepPair = '[5:44]' ;

Input(251).cellname = '100711Bc2' ; % on-off ds pair alt-v moving bar, Fred, 1+/1- yfp "okay"
Input(251).ITC18flag = 1 ;
Input(251).triggerMiss = 1 ;
Input(251).CADS = '[35:39]' ; % cell 1
Input(251).CADS = '[40:44]' ; % cell 2?
Input(251).AltVdsPair = '[55:114]' ; 

Input(252).cellname = '100711Bc2b' ; % pair above v-step data
Input(252).ITC18flag = 1 ;
Input(252).VstepPair = '[0:39]' ;

Input(253).cellname = '102111Bc1' ; % on-off ds pair alt-v moving bar, Fred, "decent" 
Input(253).ITC18flag = 1 ;
Input(253).triggerMiss = 1 ;
Input(253).CADS = '[3:7]' ; % cell 1 + yfp
Input(253).CADS = '[0:2]' ; % cell 2 - yfp
Input(253).AltVdsPair = '[13:63]' ; 

Input(254).cellname = '102111Bc1b' ; % pair above - no evidence for coupling
Input(254).ITC18flag = 1 ;
Input(254).VstepPair = '[1:170]' ; 

Input(255).cellname = '102111Bc2' ; % on-off ds pair alt-v moving bar, Fred, "not great pair - bad access on one cell though the other nice" 
Input(255).ITC18flag = 1 ;
Input(255).triggerMiss = 1 ;
Input(255).CADS = '[0:2]' ; % cell 2 - yfp
Input(255).CADS = '[3:5]' ; % cell 1 + yfp
Input(255).AltVdsPair = '[14:19,20:25,27:41,42:50]' ; 

% single cell data

Input(256).cellname = '100711Bc1' ; % on-off ds pair alt-v moving bar, Fred, 1 +yfp 1-yfp, "nice ds pair gain of second cell 20 for constant hold at end- may be messed up in file (first cell gain 10)
Input(256).ITC18flag = 1 ;
Input(256).triggerMiss = 1 ;
Input(256).CADS = '[0:3]' ; % cell 1
Input(256).AltDS = '[12:86]' ; 
Input(256).channel = 0 ;

Input(257).cellname = '100711Bc1' ; % on-off ds pair alt-v moving bar, Fred, 1 +yfp 1-yfp, "nice ds pair gain of second cell 20 for constant hold at end- may be messed up in file (first cell gain 10)
Input(257).ITC18flag = 1 ;
Input(257).triggerMiss = 1 ;
Input(257).CADS = '[0:3]' ; % cell 1
Input(257).AltDS = '[12:86]' ; 
Input(257).channel = 1 ;

Input(258).cellname = '102111Bc1' ; % on-off ds pair alt-v moving bar, Fred, "decent" 
Input(258).ITC18flag = 1 ;
Input(258).triggerMiss = 1 ;
Input(258).CADS = '[3:7]' ; % cell 1 + yfp
Input(258).CADS = '[0:2]' ; % cell 2 - yfp
Input(258).AltDS = '[13:63]' ; 
Input(258).channel = 0 ;

Input(259).cellname = '102111Bc1' ; % on-off ds pair alt-v moving bar, Fred, "decent" 
Input(259).ITC18flag = 1 ;
Input(259).triggerMiss = 1 ;
Input(259).CADS = '[3:7]' ; % cell 1 + yfp
Input(259).CADS = '[0:2]' ; % cell 2 - yfp
Input(259).AltDS = '[13:63]' ; 
Input(259).channel = 1 ;

Input(260).cellname = '102111Bc2' ; % on-off ds pair alt-v moving bar, Fred, "not great pair - bad access on one cell though the other nice" 
Input(260).ITC18flag = 1 ;
Input(260).triggerMiss = 1 ;
Input(260).CADS = '[0:2]' ; % cell 2 - yfp
Input(260).CADS = '[3:5]' ; % cell 1 + yfp
Input(260).AltDS = '[14:19,20:25,27:41,42:50]' ;
Input(260).channel = 1 ;


% 042811Fc1 on-off ds ca steps +/- gabazine 
% 042611Fc1 off s checkerboard exc and inh +/- apb/ly
% 031309E on-off +/- curare, apv and nbqx
% 031709Ec1 on-off w/apb and nbqx did not block off response

% 092309Bc1-2 on-off ds iclamp looking for spikelets with hyperpol, ttx puff and bath app

% 042010Bc1 offt, ca and wc, steps and gradings, +/-ly/apb

% 040610Bc3 on parasol, ca, gradings
% 040610Bc4 on parasol, ca and wc, gradings
% 040610Bc5 on parasol, ca, gradings
% 040610Bc6 off parasol, ca +/- ly/apb

% 060410Fc1 ds cell, wc, cs internal gaba puff and moving bars
% 060910Fc1/2 ds cell, pp gaba puff and and light stim
% 072310Fc1-3 ds cell, pp gaba puff and light steps



% Parameters assuming sample rate at 10 kHz
Parameters.PrePnts = 10000 ;    % number of points before stimuli starts
Parameters.StmPnts = 60000 ;    % number of points during stimuli
Parameters.PostPnts = 1000 ;    % number of points after stimuli ends
Parameters.STAPnts = 3000 ;     % number of points used for prespike wave forms
Parameters.DecimatePnts = 10 ;  % number of points averaged together in order to downsample prespike waveforms
Parameters.SmoothPnts = 100 ;   % number of points used to smooth spike train to make PSTH
Parameters.QuietPnts = 200 ;    % number of points of no spikes that identifiies a burst
Parameters.MinSpikeNumber = 2 ; % minimum number of spikes defining a burst
Parameters.isicut = 20 ;        % time to cut isi in half
Parameters.DClight = 1 ;        % 1=1st set of g, 2=2nd set of g
Parameters.ChopPnts = 6000 ;    % number of points used to find lin filter

Parameters.residualOption = 1 ; % 0 = calculate mean from all trials, 1 = mean of surrounding 2 trials, 2 = subtract 1 trial from nearest neighbor

Parameters.Einh = -80 ;         % inh reversal potential (mV)
Parameters.Eexc = 0 ;           % exc reversal potential (mV)
Parameters.Eleak = -70 ;        % leak g reversal potential (mV)
Parameters.cap = 0.04 ;         % cell capacitance (nF)
Parameters.Vrest = -70 ;        % initial voltage of cell (mV)      
Parameters.Vthresh = -45 ;       % spike threshold (mV)
Parameters.AbsRef = 0.002 ;     % absolute refractory period (sec) 
Parameters.RelRefTau = 0.008 ;  % relative refractory decay tau (sec)
Parameters.RelRefAmp = 4 ;      % amplitude of relative ref period (mV)
Parameters.Gleak = 1 ;          % amplitude of leak conductance (nS)

%in seconds
Parameters.IVtimes = [] ;    % points you want an IV plot

% prep matrix if it doesn't exist already
if ~exist('ForIgor')
    ForIgor=struct() ;
end

%% assess recoding quality (Cell attached only) 
if perform.CellCheckCA == 1 ;
    id = 'CA' ;
    Temp = CellCheck(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

%% assess CA burst data
if perform.BStatCA == 1 ;
    id = 'CA' ;
    Temp = BStat(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end
    
%% perfrom isi dis analysis
if perform.ISIdisCA == 1 ;
    id = 'CA' ;
    Temp = ISIdisAnalysis(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

%% perform linear nonlinear model from light to spikes
if perform.LNfiltersCA == 1 ;
    id = 'CA'
    Temp = LNfilters(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

%% perfrom variance assessment
if perform.VarGcheckExc == 1 ;
    id = 'Exc' ;
    Temp = GNoisevariance(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end
    
if perform.VarGcheckInh == 1 ;
    id = 'Inh' ;
    Temp = GNoisevariance(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

%% perform AltV analysis
if perform.AltVanalysis == 1 ;
    id = 'Alt' ;
    Temp = ValtAnalyzer(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.AltVanalysisNegCorr == 1 ;
    id = 'AltNegCorr' ;
    Temp = ValtAnalyzer(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.AltVanalysisPosCorr == 1 ;
    id = 'AltPosCorr' ;
    Temp = ValtAnalyzer(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.AltVanalysisDS == 1 ;
    id = 'AltDS' ;
    id2 = 'CADS' ;
    Temp = ValtAnalyzerDS(Input,Parameters,id,id2,A) ;
    ForIgor = mergeStruct(ForIgor, Temp) ;
end

if perform.AltInhVanalysis == 1 ;
    id = 'AltInh' ;
    Temp = ValtAnalyzer(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.AltInhVanalysisNegCorr == 1 ;
    id = 'AltNegCorrInh' ;
    Temp = ValtAnalyzer(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.AltInhVanalysisPosCorr == 1 ;
    id = 'AltPosCorrInh' ;
    Temp = ValtAnalyzer(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.AltVanalysis2 == 1 ;
    id = 'AltInh' ;
    id2 = 'Alt' ;
    Temp = ValtAnalyzer2(Input,Parameters,id,id2,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.AltVanalysisNegCorr2 == 1 ;
    id = 'AltNegCorr' ;
    Temp = ValtAnalyzer2(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.AltVanalysisPosCorr2 == 1 ;
    id = 'AltPosCorr' ;
    Temp = ValtAnalyzer2(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.AltVanalysisDS3 == 1 ;
    id = 'AltDS' ;
    Temp = ValtAnalyzerDS3(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor, Temp) ;
end

%% perform IV plot

if perform.IVSteps == 1 ;
    id = 'Steps' ;
    Temp = IVsteps(Input,Parameters,id, A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.IVStepsAPV == 1 ;
    id = 'StepsAPV' ;
    Temp = IVsteps(Input,Parameters,id, A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end 

%% perform integrate voltage model (LIF without the Gleak or Spikes) 

if perform.LIFforSimG == 1
    id = 'Alt' ;
    Temp = LIFforSimG(Input,Parameters,id, A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end 


%% perform analysis on midget or on-off ds cells dynamic clamp with alt-v g

if perform.DCmidgetAnalyzer == 1
    id = 'dcSimShuffle' ;
    Temp = DCmidgetAnalyzer(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.DCmidgetAnalyzer2 == 1
    id = 'dcSimShuffle' ;
    id2 = 'dcSimRepeat' ;
    Temp = DCmidgetAnalyzer2(Input,Parameters,id,id2,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.DCdsAnalyzer == 1 ;
    id = 'dcSimShuffle' ;
    id2 = 'dcSimRepeat' ;
    Temp = DCdsAnalyzer(Input,Parameters,id,id2,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
    % pop data
    %ForIgor.stdAll = [ForIgor.stdAll,ForIgor.(['stdAllcell',num2str(A)])] ;
    %ForIgor.ccPeakAll = [ForIgor.ccPeakAll,ForIgor.(['ccPeakAll',num2str(A)])] ;
end    

if perform.DCdsAnalyzer2 == 1 ;
    id = 'dcSimShuffle' ;
    id2 = 'dcSimRepeat' ;
    Temp = DCdsAnalyzer2(Input,Parameters,id,id2,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end    

%% perform analysis on alternating voltage exp to get time to isolate rev potential
if perform.ValtTauAlt == 1 ;
    id = 'Alt' ;
    Temp = ValtTau(Input,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.ValtTauAltDS == 1 ;
    id = 'AltDS' ;
    Temp = ValtTau(Input,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.ValtTauAltPosCorr == 1 ;
    id = 'AltPosCorr' ;
    Temp = ValtTau(Input,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.ValtTauAltNegCorr == 1 ;
    id = 'AltNegCorr' ;
    Temp = ValtTau(Input,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.ValtTauAltInh == 1 ;
    id = 'AltInh' ;
    Temp = ValtTau(Input,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

%% comparing control current step response with + drug
if perform.icomparison == 1 ;
    Temp = Icomparison(Input,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

%% ds cell attached pair data
if perform.pairCADS == 1 ;
    id = 'CADS' ;
    Temp = pairCADS(Input,Parameters,id,A) ; 
    ForIgor = mergeStruct(ForIgor,Temp) ;
end   

%% ds whole cell pair data
if perform.pairWCDSexcinh == 1 ;
end
    
%% ds alternating voltage 
if perform.pairAltDS2 == 1 ;
    id = 'AltVdsPair' ;
    Temp = pairAltDS2(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end   

%% first order LN model from light to conductance
if perform.FirstOrderLNlight2G_Exc == 1 ;
    id = 'Exc' ;
    Temp = FirstOrderLNlight2G(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.FirstOrderLNlight2G_Inh == 1 ;
    id = 'Inh' ;
    Temp = FirstOrderLNlight2G(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.FirstOrderLNlight2G_InhApb == 1 ;
    id = 'InhApb' ;
    Temp = FirstOrderLNlight2G(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

%% dynamic clamp analysis of paired ds input +/-converging +/-pairwise correlations
if perform.DCdsPairAnalyzer_141 == 1 ;
    id = 'dc141' ;
    Temp = DCdsPairAnalyzer(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
%     % pop data
%     for a = 1:3 ;
%         identifier1 = ['snCorrs',id,'Bar',num2str(a),'All'] ;
%         identifier2 = ['snCorrs',id,'Bar',num2str(a),'cell',num2str(A)] ;
%         if ~isfield(ForIgor,identifier1) ;
%             ForIgor.(identifier1) = [] ;
%         end
%         ForIgor.(identifier1) = [ForIgor.(identifier1);ForIgor.(identifier2)] ;
%     
%         identifier1 = ['stCorrs',id,'Bar',num2str(a),'All'] ;
%         identifier2 = ['stCorrs',id,'Bar',num2str(a),'cell',num2str(A)] ;
%         if ~isfield(ForIgor,identifier1) ;
%             ForIgor.(identifier1) = [] ;
%         end
%         ForIgor.(identifier1) = [ForIgor.(identifier1);ForIgor.(identifier2)] ; 
%     end
end

if perform.DCdsPairAnalyzer_157 == 1 ;
    id = 'dc157' ;
    Temp = DCdsPairAnalyzer(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
%     % pop data
%     for a = 1:3 ;
%         identifier1 = ['snCorrs',id,'Bar',num2str(a),'All'] ;
%         identifier2 = ['snCorrs',id,'Bar',num2str(a),'cell',num2str(A)] ;
%         if ~isfield(ForIgor,identifier1) ;
%             ForIgor.(identifier1) = [] ;
%         end
%         ForIgor.(identifier1) = [ForIgor.(identifier1);ForIgor.(identifier2)] ;
%     
%         identifier1 = ['stCorrs',id,'Bar',num2str(a),'All'] ;
%         identifier2 = ['stCorrs',id,'Bar',num2str(a),'cell',num2str(A)] ;
%         if ~isfield(ForIgor,identifier1) ;
%             ForIgor.(identifier1) = [] ;
%         end
%         ForIgor.(identifier1) = [ForIgor.(identifier1);ForIgor.(identifier2)] ; 
%     end
    
end

if perform.DCdsPairAnalyzer_139 == 1 ;
    id = 'dc139' ;
    Temp = DCdsPairAnalyzer(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
     % pop data
%     for a = 1:3 ;
%         identifier1 = ['snCorrs',id,'Bar',num2str(a),'All'] ;
%         identifier2 = ['snCorrs',id,'Bar',num2str(a),'cell',num2str(A)] ;
%         if ~isfield(ForIgor,identifier1) ;
%             ForIgor.(identifier1) = [] ;
%         end
%         ForIgor.(identifier1) = [ForIgor.(identifier1);ForIgor.(identifier2)] ;
%     
%         identifier1 = ['stCorrs',id,'Bar',num2str(a),'All'] ;
%         identifier2 = ['stCorrs',id,'Bar',num2str(a),'cell',num2str(A)] ;
%         if ~isfield(ForIgor,identifier1) ;
%             ForIgor.(identifier1) = [] ;
%         end
%         ForIgor.(identifier1) = [ForIgor.(identifier1);ForIgor.(identifier2)] ; 
%     end
    
end

if perform.DCdsPairRepeatAnalyzer == 1 ;
    id = 'dcRepeat' ;
    Temp = DCdsPairRepeatAnalyzer(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.DCdsPairAnalyzer3_141 == 1 ;
    id = 'dc141' ;
    Temp = DCdsPairAnalyzer3(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.DCdsPairAnalyzer3_157 == 1 ;
    id = 'dc157' ;
    Temp = DCdsPairAnalyzer3(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.DCdsPairAnalyzer3_139 == 1 ;
    id = 'dc139' ;
    Temp = DCdsPairAnalyzer3(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

%% chapter 3 of thesis analysis

if perform.chapter3ParasolAnalysis == 1 ;
    Temp = chapter3ParasolAnalysis(Input,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.chapter3ParasolDCAnalysis == 1 ;
    Temp = chapter3ParasolDCAnalysis(Input,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.chapter3MidgetAnalysis == 1 ;
    id = 'Alt' ;
    Temp = chapter3MidgetAnalysis(Input,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

end % Cell (A) loop

% % across "cell" statistics
ForIgor = PopDataConcatAndStat(ForIgor,'average') ; % get the average population data
% ForIgor = PopDataConcatAndStat(ForIgor,'sem') ; % get the standard error of the mean for population data
% % 
% %ForIgor = PopDataConcatAndStat(ForIgor,'FracChange') ; % get the fraction change for spikeCorr-ei/spikeCorr+ei
% 
% 

if perform.DCdsPairAnalyzer3_141 == 1 ;
    id = 'dc141' ;
    ForIgor = linearDiscriminantAnalysis(ForIgor,id) ;
end
  
if perform.DCdsPairAnalyzer3_139 == 1 ;
    id = 'dc139' ;
    ForIgor = linearDiscriminantAnalysis(ForIgor,id) ;
end

if perform.DCdsPairAnalyzer3_157 == 1 ;
    id = 'dc157' ;
    ForIgor = linearDiscriminantAnalysis(ForIgor,id) ;
end


% check that field names are not too long for Igor (must be less than 32 characters)
ForIgorFields = fieldnames(ForIgor) ;
for a=1:length(ForIgorFields) ; % for every field
    fieldLength(a) = length(ForIgorFields{a}) ;
end

if max(fieldLength)<32 ;
    cd ~/data_analysis/TempToIgor
    exportStructToHDF5(ForIgor,'ValtAnalyzerDS3.h5','/')
else disp('name too long for igor')
end








%% this script will analyze data from cells recorded in dynamic clamp,
% current clamp, and cell attached modes to assess basic measures of their spike output.

% JC 8/6/08


% Get the following outputs (1= get output, 0= don't) 
% basic spike and/or voltage characteristics to assess quality of recorded epochs
perform.CellCheckCA = 0 ;
perform.CellCheckIC = 0 ;
perform.CellCheckDC = 0 ;
perform.CellCheckDCnoInh = 0 ;
perform.CellCheckDCApb = 0 ;
perform.CellCheckDCnoExc = 0 ;
perform.CellCheckDCNullApb = 0 ;
perform.CellCheckDCnoExcInh = 0 ;

% interspike interval distribution
perform.ISIdisCA = 0 ;
perform.ISIdisIC = 0 ;
perform.ISIdisDC = 0 ;
perform.ISIdisDCnoInh = 0 ;
perform.ISIdisDCApb = 0 ;
perform.ISIdisDCnoExc = 0 ;
perform.ISIdisDCNullApb = 0 ;
perform.ISIdisDCnoExcInh = 0 ;

% burst analysis
perform.BStatCA = 0 ;
perform.BStatIC = 0 ;
perform.BStatDC = 1 ; 
perform.BStatDCnoInh = 1 ; 
perform.BStatDCnoExc = 0 ;
perform.BStatDCApb = 0 ;

% nearest neighbor analysis
perform.NNanalysisDC = 0 ;
perform.NNanalysisDCnoInh = 0 ;
perform.NNanalysisDCnoExc = 0 ;
perform.NNanalysisDCApb = 0 ;

% injected current and integrate and fire spike prediction
perform.LIFDC = 0 ;
perform.LIFDCnoInh = 0 ;

% plotting routines
perform.plotterDC = 0 ;
perform.plotForIgor = 0 ; % currently used as script

% using the following cells...
for A = [15] ;
%% Cells for analysis

Input(1).cellname = '070908Bc1' ; % good patch, but identity was questioned
Input(1).IC = '[88:102]' ; % only the first 5 light seeds
%Input(1).IC = '[229:238]' ; % shorter pre and post
Input(1).DC = '[251:265,305:309,315:319]' ; % with synaptic blockers
%Input(1).DC = '[104:143,169:183,214:228]' ; % without blockers
Input(1).DCnoInh = '[266:270,275:279]' ;  % with blockers
%Input(1).DCnoInh = '[144:153,184:193]' ;  % without blockers
Input(1).DCnoExc = '[285:294]' ; % with blockers
%Input(1).DCnoExc = '[159:168,194:203]' ;  % without blockers
Input(1).DCnoExcInh = '[295:304]' ;   % with blockers
%Input(1).DCnoExcInh = '[154:158,204:213]' ;   % without blockers
Input(1).DCstm = '[88:92]' ;    % stimulus used when geting G injected in DC

Input(2).cellname = '061908Bc1b' ;  % good patch but dynamic clamp was not optimized yet 
Input(2).IC = '[34:43,110:124]' ; % only first 5 light seeds
Input(2).DC = '[175:179,195:199,220:224,240:244]' ; % with synaptic blockers
%Input(2).DC = '[50:54,75:79.100:104,135:139,150:154,165:169]' ; % without blockers
Input(2).DCnoInh = '[180:184,205:209,225:229]' ;  % with blockers
%Input(2).DCnoInh = '[55:59,80:84,140:149]' ;  % without blockers
Input(2).DCnoExc = '[185:189,210:214,230:234]' ; % with blockers
%Input(2).DCnoExc = '[60:64,85:89,125:134]' ;  % without blockers
Input(2).DCnoExcInh = '[190:194]' ;   % with blockers
%Input(2).DCnoExcInh = '[65:74,90:94,155:164]' ;   % without blockers
Input(2).DCstm = '[34:38]' ;    % stimulus used when geting G injected in DC

Input(3).cellname = '061908Bc1' ; % same cell as above 
Input(3).CA = '[267:281]' ;

Input(4).cellname = '092408Bc1' ; % good patch no light response from naked retina
Input(4).DC = '[32:36,52:61,67:71,82:86,102:106]' ; % sr=10kHz unless otherwise noted
%Input(4).DC = '[122:126]' ; % sample rate = 20kHz
Input(4).DCnoInh = '[37:51,62:66,87:91,107:111]' ;
%Input(4).DCnoInh = '[127:131]' ; % sr=20kHz
Input(4).DCnoExc = '[72:76,92:96,112:116]' ;
%Input(4).DCnoExc = '[132:136]' ; %sr=20kHz
Input(4).DCnoExcInh = '[77:81,97:101,117:121]' ;
%Input(4).DCnoExcInh = '[137:141]' ; %st=20kHz

Input(5).cellname = '092408Bc2' ; % patch not great no light response from naked retina
Input(5).DC = '[20:24,35:39,50:54,60:64]' ; % amount injected current decreased between 20:39 and 50:end
Input(5).DCnoInh = '[25:29,40:44,55:59,65:69]' ;
Input(5).DCnoExc = '[30:34,75:79,85:89]' ;
Input(5).DCnoExcInh = '[45:49,80:84,90:94,105:109]' ;

Input(6).cellname = '092408Bc3' ; % good patch no light response from naked retina
Input(6).DC = '[12:16,27:36,47:51,77:86,112:121,152:156]' ; % hyperpol current injected after epoch 66, and out of blockers after epoch 71
Input(6).DCnoInh = '[17:26,37:46,52:56,87:96,122:131]' ;
Input(6).DCnoExc = '[57:61,97:106,132:141]' ;
Input(6).DCnoExcInh = '[62:66,67:71,72:76,107:111,142:151]' ;

% the cells below are using a different set of conductances than the set above
Input(7).cellname = '102108Bc1' ; % not great access (cell image kinda sparse - inj current did not adjust rest for all traces)
Input(7).DC = '[80:84]' ;
%Input(7).DC = '[118:122,133:137,143:147]' ; % inj current to adjust rest to -45 w/DC also different amp G 
Input(7).DCnoInh = '[85:89]' ;
%Input(7).DCnoInh = '[123:127,138:142]' ; % inj current to adjust rest to -45 w/DC also different amp G
Input(7).DCApb = '[90:94]' ;
%Input(7).DCApb = '[128:132,148:152]' ; % inj current to adjust rest to -45 w/ DC also different amp G
Input(7).DCNullApb = '[95:99]' ;
%Input(7).DCNullApb = '[153:157]' ; % inj current to adjust rest to -45 w/ DC also different amp G 
Input(7).DCnoExc = '[100:104]' ;
%Input(7).DCnoExc = '[105:108]' ; % inj current to adjust rest to -45 w/ DC also different amp G

Input(8).cellname = '102108Bc2' ; % decent patch
Input(8).DC = '[78:82,108:112,128:132,153:157]' ;
Input(8).DCnoExc = '[83:87,148:152]' ;
Input(8).DCnoInh = '[88:92,113:117,133:137,158:162]' ;
Input(8).DCApb = '[93:97,118:122,138:142,163:167]' ;
Input(8).DCNullApb = '[98:102,123:127]' ;

Input(9).cellname = '102108Bc3' ; % decent patch (past epoch 160 out of blockers)
Input(9).DC = '[63:72,93:97,113:117,150:154,175:179]' ;
Input(9).DCnoInh = '[73:77,98:102,123:127,155:159,180:184]' ;
Input(9).DCApb = '[78:82,103:107,160:164,185:189]' ;
Input(9).DCNullApb = '[83:87,108:112,165:169,190:194]' ;
Input(9).DCnoExc = '[118:122,170:174]' ;

Input(10).cellname = '102108Bc4' ; % decent patch (more epochs but ossillations present)
Input(10).DC = '[31:35,56:60]' ;
Input(10).DCnoInh = '[36:40,61:65]' ;
Input(10).DCApb = '[41:45]' ;
Input(10).DCNullApb = '[46:50]' ;
Input(10).DCnoExc = '[51:55]' ;

Input(11).cellname = '102108Bc5' ; % not a good patch
Input(11).DC = '[24:28]' ;
Input(11).DCnoInh = '[29:33]' ;
Input(11).DCApb = '[34:38]' ;

%Input(12).cellname = '102108Bc6' ; % bad access

Input(12).cellname = '102208Bc1' ; % decent patch but has very low firing rate (non labeled are 20/30 G amp)
Input(12).DC = '[60:64]' ; % 30/45 G amp
%Input(12).DC = '[96:100,121:125,146:150,171:175]' ;
Input(12).DCnoInh = '[65:69]' ; % 30/45 G amp
%Input(12).DCnoInh = '[101:105,126:130,151:155,176:180]' ;
Input(12).DCApb = '[76:80]' ; % 30/45 G amp
%Input(12).DCApb = '[106:110,136:140,156:160,181:185]' ;
Input(12).DCNullApb = '[81:85]' ; % 30/45 G amp
%Input(12).DCNullApb = '[111:115,131:135,161:165]' ;
Input(12).DCnoExc = '[86:90]' ; % 30/45 G amp
%Input(12).DCnoExc = '[116:120,141:145,166:170]' ; 

Input(13).cellname = '102208Bc2' ; % decent patch but cell seemed to have depolarized rest
Input(13).DC = '[45:49,70:74,95:99]' ;
Input(13).DCnoInh = '[50:54,75:79]' ;
Input(13).DCApb = '[55:59,80:84]' ;
Input(13).DCNullApb = '[60:64,85:89]' ;
Input(13).DCnoExc = '[65:69,90:94]' ;

Input(14).cellname = '100808Ac1' ; %freds dc recording in I-clamp fast
Input(14).DC = '[0:14]' ;

Input(15).cellname = '011409Fc4' ; % good patch used I-clamp fast
%Input(15).DC = '[46:50,52:56]' ; % 10/15 (everything else 20/20) 
Input(15).DC = '[57:61,77:81,112:116,137:141]' ;
Input(15).DCnoInh = '[62:66,82:86,117:121]' ;
Input(15).DCnoExc = '[67:71,87:91,122:126]' ;
Input(15).DCApb = '[72:76,92:96,127:131]' ;
Input(15).DCnoExcInh = '[97:101]' ;
Input(15).DCNullApb = '[102:106,107:111,132:136]' ;

Input(16).cellname = '011409Fc5' ; % used I-clamp fast but ocillations seem present
%Input(16).DC = '[210:214]' ; % 40/60 
%Input(16).DC = '[215:219]' ; %10/15
%Input(16).DC = '[220:224]' ; % 20/30
Input(16).DC = '[100:104,125:129,150:154,175:179,205:209]' ; %30/45 all below also
Input(16).DCnoInh = '[105:109,130:134,155:159,180:184]' ;
Input(16).DCnoExc = '[110:114,140:144,165:169,195:199]' ;
Input(16).DCApb = '[115:119,185:189]' ;
Input(16).DCnoExcInh = '[145:149,170:174,200:204]' ;
Input(16).DCNullApb = '[120:124,135:139,160:164,190:194]' ;

Input(17).cellname = '012809Bc3' ; % Petris dc cell filtering
Input(17).DC = '[106:135]';
Input(17).DCnoInh = '[]';  %this has inhibition but it a different filter
Input(17).DCnoExc = '[]'; %this has exc but its a different filter

% Parameters assuming sample rate at 10 kHz
Parameters.PrePnts = 10000 ;    % number of points before stimuli starts
Parameters.StmPnts = 60000 ;    % number of points during stimuli
Parameters.PostPnts = 1000 ;    % number of points after stimuli ends
Parameters.STAPnts = 3000 ;     % number of points used for prespike wave forms
Parameters.DecimatePnts = 10 ;  % number of points averaged together in order to downsample prespike waveforms
Parameters.SmoothPnts = 100 ;   % number of points used to smooth spike train to make PSTH
Parameters.QuietPnts = 200 ;    % number of points used to identify quiet time before burst
Parameters.MinSpikeNumber = 2 ; % minimum number of spikes in a burst
Parameters.WCSpikeThresh = -10 ; % spike threshold to detect whole cell spikes
Parameters.isicut = 10 ; %ms    % where the isi should be cut in half

% prep matrix if it doesn't exist already
if ~exist('ForIgor')
    ForIgor=struct() ;
end

%% assess recoding quality 
if perform.CellCheckCA == 1 ;
    id = 'CA' ;
    Temp = CellCheck(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.CellCheckIC == 1 ;
    id = 'IC' ;
    Temp = CellCheck(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.CellCheckDC == 1 ;
    id = 'DC' ;
    Temp = CellCheck(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.CellCheckDCnoInh == 1 ;
    id = 'DCnoInh' ;
    Temp = CellCheck(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.CellCheckDCApb == 1 ;
    id = 'DCApb' ;
    Temp = CellCheck(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.CellCheckDCnoExc == 1 ;
    id = 'DCnoExc' ;
    Temp = CellCheck(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.CellCheckDCNullApb == 1 ;
    id = 'DCNullApb' ;
    Temp = CellCheck(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.CellCheckDCnoExcInh == 1 ;
    id = 'DCnoExcInh' ;
    Temp = CellCheck(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end


%% perfrom isi dis analysis
if perform.ISIdisCA == 1 ;
    id = 'CA' ;
    Temp = ISIdisAnalysis(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.ISIdisIC == 1 ;
    id = 'IC' ;
    Temp = ISIdisAnalysis(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.ISIdisDC == 1 ;
    id = 'DC' ;
    Temp = ISIdisAnalysis(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.ISIdisDCnoInh == 1 ;
    id = 'DCnoInh' ;
    Temp = ISIdisAnalysis(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.ISIdisDCApb == 1 ;
    id = 'DCApb' ;
    Temp = ISIdisAnalysis(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.ISIdisDCnoExc == 1 ;
    id = 'DCnoExc' ;
    Temp = ISIdisAnalysis(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.ISIdisDCNullApb == 1 ;
    id = 'DCNullApb' ;
    Temp = ISIdisAnalysis(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.ISIdisDCnoExcInh == 1 ;
    id = 'DCnoExcInh' ;
    Temp = ISIdisAnalysis(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

%%
if perform.BStatCA == 1 ;    
    id ='CA' ;
    Temp = BStat(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end
    
if perform.BStatIC == 1 ;
    id ='IC' ;
    Temp = BStat(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.BStatDC == 1 ;
    id ='DC' ;
    Temp = BStat(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.BStatDCnoInh == 1 ;
    id ='DCnoInh' ;
    Temp = BStat(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.BStatDCnoExc == 1 ;
    id ='DCnoExc' ;
    Temp = BStat(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.BStatDCApb == 1 ;
    id ='DCApb' ;
    Temp = BStat(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

%% Nearest neighbor analysis

if perform.NNanalysisDC == 1 ;
    id = 'DC' ;
    Temp = NNanalysis(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.NNanalysisDCnoInh == 1 ;
    id = 'DCnoInh' ;
    Temp = NNanalysis(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.NNanalysisDCnoExc == 1 ;
    id = 'DCnoExc' ;
    Temp = NNanalysis(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.NNanalysisDCApb == 1 ;
    id = 'DCApb' ;
    Temp = NNanalysis(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

%% injected current and integrate and fire spike prediction

if perform.LIFDC == 1 ;
    id = 'DC' ;
    Temp = dcLIF(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.LIFDCnoInh == 1 ;
    id = 'DCnoInh' ;
    Temp = dcLIF(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end


%% plotting routines

if perform.plotterDC == 1 ;
    id = 'DC' ;
    ploted = plotterDC(Input,Parameters,id,A) ; % plots DC V data and normalized G injected
end

end % end cell for loop

if perform.plotForIgor == 1;
    ploted = plotForIgor(ForIgor) ; % this is a collection of ploting routines which will plot data in the ForIgor matrix 
end
    

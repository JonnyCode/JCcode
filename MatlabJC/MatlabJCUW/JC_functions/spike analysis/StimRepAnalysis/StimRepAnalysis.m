%% this script will analyze data from cells recorded in dynamic clamp,
% current clamp, and cell attached modes to assess how the spikes represent the stimulus.

% JC 8/6/08


% Get the following outputs (1= get output, 0= don't) 
% Prespike stimulus (STA), linear and nonlinear filters (LN model without
% spike generation)
perform.LNCA = 0 ;         % cell attached recording
perform.LNIC = 0 ;         % Current clamp
perform.LNDC = 0 ;         % dynamic clamp wiht both exc and inh
perform.LNDCnoInh = 0 ;    % dynamic clamp with exc but no timevarying inh  
perform.LNDCnoExc = 0 ;    % dynamic clamp with inh but no timevarying exc 
perform.LNDCnoExcInh = 0 ; % dynamic clamp with neither timevarying exc or inh

% prespike stim +/- inh
perform.NewSpkSTA = 0 ;

% Spike representation of stim (reverse LN model)
perform.StmRepDC = 0 ;
perform.StmRepDCnoInh = 0 ;
perform.StmRepDCnoExc = 0 ;


% prep matrix if it doesn't exist already
if ~exist('ForIgor')
    ForIgor=struct() ;
end

% using the following cells...
for A = [4] ;
%% Cells for analysis

Input(1).cellname = '070908Bc1' ; % good patch, but identity was questioned
Input(1).IC = '[88:102,229:238]' ;
Input(1).DC = '[251:265,305:309,315:329]' ; % with synaptic blockers
%Input(1).DC = '[104:143,169:183,214:228]' ; % without blockers
Input(1).DCnoInh = '[266:270,275:279]' ;  % with blockers
%Input(1).DCnoInh = '[144:153,184:193]' ;  % without blockers
Input(1).DCnoExc = '[285:294]' ; % with blockers
%Input(1).DCnoExc = '[159:168,194:203]' ;  % without blockers
Input(1).DCnoExcInh = '[295:304]' ;   % with blockers
%Input(1).DCnoExcInh = '[154:158,204:213]' ;   % without blockers
Input(1).DCstm = '[88:92]' ;    % stimulus used when geting G injected in DC

Input(2).cellname = '061908Bc1b' ;  % good patch but dynamic clamp was not optimized yet 
Input(2).IC = '[34:43,110:124]' ;
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
Input(4).DC = '[32:36,52:61,67:71,82:86,102:106,149:168]' ; % sr=10kHz unless otherwise noted
%Input(4).DC = '[122:126]' ; % sample rate = 20kHz
Input(4).DCnoInh = '[37:51,62:66,87:91,107:111,169:173]' ;
%Input(4).DCnoInh = '[127:131]' ; % sr=20kHz
Input(4).DCnoExc = '[72:76,92:96,112:116]' ;
%Input(4).DCnoExc = '[132:136]' ; %sr=20kHz
Input(4).DCnoExcInh = '[77:81,97:101,117:121]' ;
%Input(4).DCnoExcInh = '[137:141]' ; %st=20kHz

Input(5).cellname = '092408Bc2' ; % patch not great no light response from naked retina
Input(5).DC = '[20:24,35:39,50:54,60:64,95:104]' ; % amount injected current decreased between 20:39 and 50:end
Input(5).DCnoInh = '[25:29,40:44,55:59,65:69]' ;
Input(5).DCnoExc = '[30:34,75:79,85:89]' ;
Input(5).DCnoExcInh = '[45:49,80:84,90:94,105:109]' ;

Input(6).cellname = '092408Bc3' ; % good patch no light response from naked retina
Input(6).DC = '[12:16,27:36,47:51,77:86,112:121,152:156]' ; % hyperpol current injected after epoch 66, and out of blockers after epoch 71
Input(6).DCnoInh = '[17:26,37:46,52:56,87:96,122:131]' ;
Input(6).DCnoExc = '[57:61,97:106,132:141]' ;
Input(6).DCnoExcInh = '[62:66,67:71,72:76,107:111,142:151]' ;

% Parameters assuming sample rate at 10 kHz
Parameters.PrePnts = 10000 ;    % number of points before stimuli starts
Parameters.StmPnts = 60000 ;    % number of points during stimuli
Parameters.PostPnts = 1000 ;    % number of points after stimuli ends
Parameters.STAPnts = 3000 ;     % number of points used for prespike wave forms
Parameters.DecimatePnts = 10 ;  % number of points averaged together in order to downsample prespike waveforms
Parameters.SmoothPnts = 100 ;   % number of points used to smooth spike train to make PSTH
Parameters.QuietPnts = 200 ;    % number of points used to identify quiet time before burst
Parameters.MinSpikeNumber = 2 ; % minimum number of spikes in a burst

%% Linear Nonlinear model

if perform.LNCA == 1 ; 
    id = 'CA' ;
    Temp = LNfilters(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.LNIC == 1 ; 
    id = 'IC' ;
    Temp = LNfilters(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.LNDC == 1 ; 
    id = 'DC' ;
    Temp = LNfilters(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.LNDCnoInh == 1 ; 
    id = 'DCnoInh' ;
    Temp = LNfilters(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.LNDCnoExc == 1 ; 
    id = 'DCnoExc' ;
    Temp = LNfilters(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.LNDCnoExcInh == 1 ; 
    id = 'DCnoInhExc' ;
    Temp = LNfilters(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

%% How well do the spikes represent the stimulus?
if perform.StmRepDC == 1 ;
    id = 'DC' ;
    Temp = LNstmRep(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.StmRepDCnoInh == 1 ;
    id = 'DCnoInh' ;
    Temp = LNstmRep(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.StmRepDCnoExc == 1 ;
    id = 'DCnoExc' ;
    Temp = LNstmRep(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

%% New vs missing spikes

if perform.NewSpkSTA == 1 ;
    id1 = 'DC' ;
    id2 = 'DCnoInh' ;
    Temp = NewSpksAnalysis(Input,Parameters,id1,id2,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

end % end for loop

% extra figures
% figure
% plot(ForIgor.DisXDC1,ForIgor.stmDisDCnoInh1/sum(ForIgor.stmDisDCnoInh1),'k')
% hold on
% plot(ForIgor.DisXDC1,ForIgor.predDisDC1/sum(ForIgor.predDisDC1),'b')
% plot(ForIgor.DisXDC1,ForIgor.predDisDCnoInh1/sum(ForIgor.predDisDCnoInh1),'r')
% plot(ForIgor.DisXDC1,ForIgor.predDisDCnoExc1/sum(ForIgor.predDisDCnoExc1),'g')
% 
% for a=[1,4:6]
% figure(a)
% subplot(2,3,1)
% plot(-[1:length(ForIgor.(['STADC',num2str(a)]))],ForIgor.(['STADC',num2str(a)]))
% hold on
% plot(-[1:length(ForIgor.(['STADCnoInh',num2str(a)]))],ForIgor.(['STADCnoInh',num2str(a)]),'r')
% axis([-9000 0 -50 150])
% xlabel('sample points')
% ylabel('Light')
% title('STA')
% 
% subplot(2,3,4)
% plot(-[1:length(ForIgor.(['STADC',num2str(a)]))],ForIgor.(['STADC',num2str(a)])/max(ForIgor.(['STADC',num2str(a)])))
% hold on
% plot(-[1:length(ForIgor.(['STADCnoInh',num2str(a)]))],ForIgor.(['STADCnoInh',num2str(a)])/max(ForIgor.(['STADCnoInh',num2str(a)])),'r')
% axis([-9000 0 -.5 1])
% title('STAnorm')
% 
% subplot(2,3,2)
% plot(-[1:length(ForIgor.(['LinFilterDC',num2str(a)]))],ForIgor.(['LinFilterDC',num2str(a)]))
% hold on
% plot(-[1:length(ForIgor.(['LinFilterDCnoInh',num2str(a)]))],ForIgor.(['LinFilterDCnoInh',num2str(a)]),'r')
% axis([-9000 0 -.000005 .00001])
% title('linear filter')
% 
% subplot(2,3,5)
% plot(-[1:length(ForIgor.(['LinFilterDC',num2str(a)]))],ForIgor.(['LinFilterDC',num2str(a)])/max(ForIgor.(['LinFilterDC',num2str(a)])))
% hold on
% plot(-[1:length(ForIgor.(['LinFilterDCnoInh',num2str(a)]))],ForIgor.(['LinFilterDCnoInh',num2str(a)])/max(ForIgor.(['LinFilterDCnoInh',num2str(a)])),'r')
% axis([-9000 0 -.5 1])
% title('linear filter norm')
% 
% subplot(2,3,3)
% plot(ForIgor.(['generatorBinsDC',num2str(a)]),ForIgor.(['SpksPerPntsDC',num2str(a)]),'b')
% hold on
% plot(ForIgor.(['generatorBinsDCnoInh',num2str(a)]),ForIgor.(['SpksPerPntsDCnoInh',num2str(a)]),'r')
% xlabel('generator bins')
% ylabel('mean spikes/pnt')
% title('nonlinear')
% end
% 
% 

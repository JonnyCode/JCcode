function ForIgor = chapter3MidgetAnalysis(Input,id,A) 

% thesis data and analysis for chapter 3

[fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/primate/',Input(A).cellname]);
 
epochs = str2num(Input(A).(id)) ;
round = 0 ;
for a = 1:3:length(epochs) ; % for each spike epoch
    round = round +1 ;

    [dataExc(round,:), error] = ITCReadEpoch(epochs(a), 0, fp) ;    % get data

    [dataInh(round,:), error] = ITCReadEpoch(epochs(a+1), 0, fp) ;    % get data
    
    [dataAltV(round,:), error] = ITCReadEpoch(epochs(a+2), 0, fp) ;    % get data
    
    [lightCommand(round,:), error] = ITCReadEpochStm(epochs(a+2), 0,fp); % get light
    [voltageCommand(round,:), error] = ITCReadEpochStm(epochs(a+2), 1,fp); % get voltage command
    
end

[SI, error] = ITCGetSamplingInterval(epochs(1), fp); % get sampling interval
SI = SI * 1e-6; % Sampling interval in sec
if Input(A).ITC18flag == 1 ;
    SI = SI*1.25 ;
end

lightCommand = lightCommand(:,1:length(dataExc)) ;
lightCommand(lightCommand<0)=0 ;

voltageCommand = voltageCommand(:,1:length(dataExc)) ;

cyclepnts = 100 ; % number of sample points in a cycle, pnts between leaving hold1 and returning (Also gets rid of first cycle)    
FirstAltPnt = (cyclepnts/2)+1 ; % first sample point you want to plot after begining of step from alternation 
LastAltPnt = FirstAltPnt ;  % last "                                                                    "

[prePnts, error] = ITCGetStmPrePts(epochs(1), 0, 0, fp) ; % points collected beyond which data is not worth analyzing
[postPnts, error] = ITCGetStmTailPts(epochs(1), 0, 0, fp) ;

samplerate = 1/SI(1) ; % Hz at which data was collected
time = [SI(1):SI(1):SI(1)*length(dataExc(1,:))] ; % time vector in seconds

% WINDOWING ALTERNATING V DATA 
Alt_Exc = NaN(size(dataAltV)) ; % make a vector of NaNs that will serve as base for ploting alternating coductances at hold1
Alt_Inh = NaN(size(dataAltV)) ; % make a vector of NaNs that will serve as base for ploting alternating coductances at hold2_

for a = FirstAltPnt:LastAltPnt ;        % for each data point per cycle we want to plot
    Alt_Exc(:,[a:cyclepnts:end]) = dataAltV(:,[a:cyclepnts:end]) ; 
    Alt_Inh(:,[a+cyclepnts/2:cyclepnts:end]) = dataAltV(:,[a+cyclepnts/2:cyclepnts:end]) ; 
end 

% interpolate alternating current data
for a = 1:size(dataExc,1) ; % for each trial
     b = find(isnan(Alt_Exc(a,:)) == 0) ;     % find all the indices that are not nans
    Alt_ExcInt(a,:) = interp1(b,Alt_Exc(a,b),[1:length(Alt_Exc)],'linear','extrap') ; % interpolate to find values that were not sampled
    clear b
    
    b = find(isnan(Alt_Inh(a,:)) == 0) ;     % find all the indices that are not nans
    Alt_InhInt(a,:) = interp1(b,Alt_Inh(a,b),[1:length(Alt_Inh)],'linear','extrap') ; % interpolate to find values that were not sampled
    clear b
end

% low pass filter to remove electrical crap
dataExc_lpf = lowPassFilter(dataExc, samplerate, 5000) ; %(signal,samplerate,cutoff frequ (hz))
dataInh_lpf = lowPassFilter(dataInh, samplerate, 5000) ; 
Alt_ExcInt_lpf = lowPassFilter(Alt_ExcInt, samplerate, 5000) ;
Alt_InhInt_lpf = lowPassFilter(Alt_InhInt, samplerate, 5000) ;

% high pass filter to remove slow drift 
% greg S. sent this to me to help implement a butterworth and avoid ringing
F=1 ; % filter cuttoff
Wn = F*SI(1); %normalized frequency cutoff
[z, p, k] = butter(1,Wn,'high'); %
[sos,g]=zp2sos(z,p,k); 
myfilt=dfilt.df2sos(sos,g);

dataExc_hpf = filter(myfilt,dataExc_lpf')'; % filter implementation
dataInh_hpf = filter(myfilt,dataInh_lpf')'; 
Alt_ExcInt_hpf = filter(myfilt,Alt_ExcInt_lpf')'; 
Alt_InhInt_hpf = filter(myfilt,Alt_InhInt_lpf')';


% change to conductances and subtract off means 
for a = 1:size(dataExc,1) ; % for each trial
    GExc_hpf(a,:) = dataExc_hpf(a,:)/-61 - mean(dataExc_hpf(a,prePnts:end-postPnts)/-61) ; % get conductance from stable currrents
    GInh_hpf(a,:) = dataInh_hpf(a,:)/61 - mean(dataInh_hpf(a,prePnts:end-postPnts)/61) ; 
    GAlt_ExcInt_hpf(a,:) = Alt_ExcInt_hpf(a,:)/-61 - mean(Alt_ExcInt_hpf(a,prePnts:end-postPnts)/-61) ; % get conductance from alt current
    GAlt_InhInt_hpf(a,:) = Alt_InhInt_hpf(a,:)/61 - mean(Alt_InhInt_hpf(a,prePnts:end-postPnts)/61) ; %#ok<*AGROW>

end

offsetExc = min(min([GExc_hpf(:,prePnts:end-postPnts),GAlt_ExcInt_hpf(:,prePnts:end-postPnts)])) ;
offsetInh = min(min([GInh_hpf(:,prePnts:end-postPnts),GAlt_InhInt_hpf(:,prePnts:end-postPnts)])) ;
% ofset G (assumes all g have same mean and 1 min) 
for a = 1:size(dataExc,1) ; % for each trial
    GExc_hpf(a,:) = GExc_hpf(a,:) - offsetExc ; % offsets
    GInh_hpf(a,:) = GInh_hpf(a,:) - offsetInh ; 
    GAlt_ExcInt_hpf(a,:) = GAlt_ExcInt_hpf(a,:) - offsetExc ; 
    GAlt_InhInt_hpf(a,:) = GAlt_InhInt_hpf(a,:) - offsetInh ;
    
end

% get mean data
GAlt_ExcInt_hpf_Mean = mean(GAlt_ExcInt_hpf) ;
GAlt_InhInt_hpf_Mean = mean(GAlt_InhInt_hpf) ;

GExc_hpf_Mean = mean(GExc_hpf) ;
GInh_hpf_Mean = mean(GInh_hpf) ;


% g prep
GExc = GExc_hpf_Mean(prePnts:end-postPnts) ;
GInh = GInh_hpf_Mean(prePnts:end-postPnts) ;

GExc = GExc - min(GExc) ;
GInh = GInh - min(GInh) ;

ccG = xcov(GExc,GInh,'coeff') ;
ccTime = ([1:length(ccG)] - (length(ccG)+1)/2)*SI ;


params.Eexc = 0 ;          % reversal potential for excitatory current
params.Einh = -70;
params.Eleak = -55 ;         % reversal potential of leak conductance (this can make cell spike spontaneously)

params.cap = 0.045 ; %nF                          % Set capacitance of cell
params.Vrest = params.Eleak ;                 % initial potential of cell (resting potential mV) 
params.Vthresh = -40 ;     % spike threshold
params.AbsRef = .002 ; %sec       % absolute refractory period 

params.RelRefTau = .01 ; %sec              % decay time constant of relative refractory period
params.RelRefAmp = 4 ; %mV                  % amplitude of relative refractory period threshold change  

Iadd = 0 ; % pA
Gleak = 1 ; % nS 
Gleak_range = [Gleak:.1:5] ;

for a=1:size(GExc_hpf,1) ;
    LIFv_Control(a,:) = LIFmodelGplusI(GExc_hpf(a,:),GInh_hpf(a,:),Gleak,Iadd,samplerate,params) ;
    LIFv_minusInh(a,:) = LIFmodelGplusI(GExc_hpf(a,:),ones(1,length(GInh_hpf))*mean(GInh_hpf(a,prePnts-1000:prePnts)),Gleak,Iadd,samplerate,params) ;

end

LIFv_Control_spikeTrain = zeros(size(LIFv_Control)) ;
LIFv_minusInh_spikeTrain = zeros(size(LIFv_minusInh)) ;

LIFv_Control_spikeTrain(LIFv_Control==50) = 1 ;
LIFv_minusInh_spikeTrain(LIFv_minusInh==50) = 1 ;

for a=1:size(LIFv_Control_spikeTrain,1) ;
    LIFv_Control_spikeTimes{a} = time(LIFv_Control(a,:)==50) ;
    LIFv_minusInh_spikeTimes{a} = time(LIFv_minusInh(a,:)==50) ;
end

% stas
staTime = .3 ; %sec
staLength = staTime/SI ; % pnts

staControl = staFinder(lightCommand(:,prePnts:end),LIFv_Control_spikeTrain(:,prePnts:end),staLength) ;
staminusInh = staFinder(lightCommand(:,prePnts:end),LIFv_minusInh_spikeTrain(:,prePnts:end),staLength) ;

staTime = [-staLength+1:0]*SI ;

% figure
% plot(time,GExc_hpf_Mean,'b--')
% hold on
% plot(time,GInh_hpf_Mean,'r--')
% plot(time,GExc_hpf(1,:))
% plot(time,GInh_hpf(1,:),'r')
% 
% 
% figure
% plot(ccTime,ccG)

% figure
% plot(time,LIFv_Control(1,:),'b')
% hold on
% plot(time,LIFv_minusInh(1,:),'r')
% 
% figure
% for a=1:size(LIFv_Control_spikeTrain,1);
%     for b=1:length(LIFv_Control_spikeTimes{a}) ;
%         plot([LIFv_Control_spikeTimes{a}(b),LIFv_Control_spikeTimes{a}(b)],[a-1,a],'b')
%         hold on
%     end
%     for b=1:length(LIFv_minusInh_spikeTimes{a}) ;
%         plot([LIFv_minusInh_spikeTimes{a}(b),LIFv_minusInh_spikeTimes{a}(b)],[length(LIFv_Control_spikeTimes)+a-1,length(LIFv_Control_spikeTimes)+a],'r')
%         hold on
%     end
% end
% 
% figure
% plot(staTime,staControl,'b') ;
% hold on
% plot(staTime,staminusInh,'r') ;


% for igor

% % example g data
% identifier = ['GexcExample',num2str(A)] ;
% ForIgor.(identifier) = GExc_hpf(1,:) ; 
% 
% identifier = ['GinhExample',num2str(A)] ;
% ForIgor.(identifier) = GInh_hpf(1,:) ; 
% 
% identifier = ['lightExample',num2str(A)] ;
% ForIgor.(identifier) = lightCommand(1,:) ; 
% 
% identifier = ['timeExample',num2str(A)] ;
% ForIgor.(identifier) = time ; 
% 
% % cross correlation
% identifier = ['cc',num2str(A)] ;
% ForIgor.(identifier) = ccG ;
% 
% identifier = ['ccTime',num2str(A)] ;
% ForIgor.(identifier) = ccTime ;

% example LIF data
identifier = ['LIFvControl',num2str(A)] ;
ForIgor.(identifier) = LIFv_Control(1,:) ; 

identifier = ['LIFvMinusInh',num2str(A)] ;
ForIgor.(identifier) = LIFv_minusInh(1,:) ; 

% example LIF spike trains
LIFv_Control_spikeTrain(1,(LIFv_Control_spikeTrain(1,:)==0)) = nan ;
LIFv_minusInh_spikeTrain(1,(LIFv_minusInh_spikeTrain(1,:)==0)) = nan ;

identifier = ['strainControlExample',num2str(A)] ;
ForIgor.(identifier) = LIFv_Control_spikeTrain(1,:) ; 

identifier = ['strainMinusInhExample',num2str(A)] ;
ForIgor.(identifier) = LIFv_minusInh_spikeTrain(1,:)*2 ; 


% sta
identifier = ['staTime',num2str(A)] ;
ForIgor.(identifier) = staTime ;

identifier = ['staControl',num2str(A)] ;
ForIgor.(identifier) = staControl ;

identifier = ['staMinusInh',num2str(A)] ;
ForIgor.(identifier) = staminusInh ;




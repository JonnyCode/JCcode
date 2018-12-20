function ForIgor = ValtAnalyzer(Input,Parameters,id,A) ; 
% this function will analyze alternate voltage exp
% JC 2/9/09 adapted from earlier version

% 11/2/09 changed highpass filter to buterworth and added ITC flag
% 11/20/09 rearranged code to filter current and change I to g conversion
% and added a bunch of ForIgors, sampled last pnt instead of 2nd to last,
% fixed resample sh
% 2/1/10 added dot product check to assess slow drift
% 8/4/10 added residual options, ac of singlehold res and cc error bars  


% Cells for analysis

% Input(1).cellname = '101907Bc1' ; % good patch, but identity was questioned
% Input(1).AltV = '[140,146,152,158,164]' ; % alternating voltages
% Input(1).HoldExc = '[142,148,154,160,166]' ; % Exc isolation hold
% Input(1).HoldInh = '[144,150,156,162,168]' ; % inh isolation hold
% 
% Input(2).cellname = '020609Ec1' ;
% Input(2).AltV = '[182:3:194]' ;
% Input(2).HoldExc = '[180:3:194]' ;
% Input(2).HoldInh = '[181:3:194]' ;

residualOption = Parameters.residualOption ;

% get data
try
    [fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/primate/',Input(A).cellname]);
catch
    [fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/mouse/',Input(A).cellname]);
end

epochs = str2num(Input(A).(id)) ;
round = 0 ;
for a = 1:3:length(epochs) ; % for each spike epoch
    round = round +1 ;

    [dataExc(round,:), error] = ITCReadEpoch(epochs(a), 0, fp) ;    % get data

    [dataInh(round,:), error] = ITCReadEpoch(epochs(a+1), 0, fp) ;    % get data
    
    [dataAltV(round,:), error] = ITCReadEpoch(epochs(a+2), 0, fp) ;    % get data
    
    [lightCommand(round,:), error] = ITCReadEpochStm(epochs(a+2), 0,fp); % get light
    [voltageCommand(round,:), error] = ITCReadEpochStm(epochs(a+2), 1,fp); % get voltage command
    
    [SI(round), error] = ITCGetSamplingInterval(epochs(a+2), fp); % get sampling interval
    SI(round) = SI(round) * 1e-6; % Sampling interval in sec
end

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

% calculate residuals by specified method (this was added 8/4/10)
if residualOption == 0 ;
    % get residuals standard calculation
    for a= 1:size(GExc_hpf,1) ;
        GExc_hpf_Res(a,:) =  GExc_hpf(a,:) - GExc_hpf_Mean ;
        GInh_hpf_Res(a,:) =  GInh_hpf(a,:) - GInh_hpf_Mean ;
        GAlt_ExcInt_hpf_Res(a,:) = GAlt_ExcInt_hpf(a,:) - GAlt_ExcInt_hpf_Mean ;
        GAlt_InhInt_hpf_Res(a,:) = GAlt_InhInt_hpf(a,:) - GAlt_InhInt_hpf_Mean ;
    end

elseif residualOption == 1 ;
    % get residual nearest neighbor calculation (added option 8/4/10)
    for a= 2:size(GExc_hpf,1)-1 ;
        GExc_hpf_Res(a-1,:) =  GExc_hpf(a,:) - (GExc_hpf(a-1,:)+GExc_hpf(a+1,:))/2 ; %ERROR FOUND 10/13/10 was written this way before: GExc_hpf(a+1,:)+GExc_hpf(a+1,:))/2
        GInh_hpf_Res(a-1,:) =  GInh_hpf(a,:) - (GInh_hpf(a-1,:)+GInh_hpf(a+1,:))/2 ;
        GAlt_ExcInt_hpf_Res(a-1,:) = GAlt_ExcInt_hpf(a,:) - (GAlt_ExcInt_hpf(a-1,:)+GAlt_ExcInt_hpf(a+1,:))/2 ;
        GAlt_InhInt_hpf_Res(a-1,:) = GAlt_InhInt_hpf(a,:) - (GAlt_InhInt_hpf(a-1,:)+GAlt_InhInt_hpf(a+1,:))/2 ;
    end

elseif residualOption == 2 ;    
    % get residual nearest neighbor calculation (added option 8/4/10) - this one uses all trials
    for a= 1:size(GExc_hpf,1) ;

        if rem(a,2)~=0 ; % if its odd
            b = a+1 ;
        else             % if its even
            b = a-1 ;
        end

        if rem(size(GExc_hpf,1),2)~=0 & size(GExc_hpf,1)-1==a; % if the number of trials is odd and its the sencond to last trial
            b = a+1 ;
        end

        if rem(size(GExc_hpf,1),2)~=0 & size(GExc_hpf,1)==a; % if the number of trials is odd and its the last trial
            b = a-2 ;
        end

        GExc_hpf_Res(a,:) =  GExc_hpf(a,:) - GExc_hpf(b,:) ;
        GInh_hpf_Res(a,:) =  GInh_hpf(a,:) - GInh_hpf(b,:);
        GAlt_ExcInt_hpf_Res(a,:) = GAlt_ExcInt_hpf(a,:) - GAlt_ExcInt_hpf(b,:) ;
        GAlt_InhInt_hpf_Res(a,:) = GAlt_InhInt_hpf(a,:) - GAlt_InhInt_hpf(b,:) ;
    end
    clear b

end    
    
% get dot product to check for slow drift
for a= 1:size(GExc_hpf,1) ;
    GExc_hpf_dot(a) =  mean(GExc_hpf(a,:).*GExc_hpf_Mean)/mean(GExc_hpf_Mean.^2) ;
    GInh_hpf_dot(a) =  mean(GInh_hpf(a,:).*GInh_hpf_Mean)/mean(GInh_hpf_Mean.^2) ;
    GAlt_ExcInt_hpf_dot(a) = mean(GAlt_ExcInt_hpf(a,:).*GAlt_ExcInt_hpf_Mean)/mean(GAlt_ExcInt_hpf_Mean.^2) ;
    GAlt_InhInt_hpf_dot(a) = mean(GAlt_InhInt_hpf(a,:).*GAlt_InhInt_hpf_Mean)/mean(GAlt_InhInt_hpf_Mean.^2) ;
end

% get cross correlations (avoid first 3sec because of oscillations from highpass filtering) 
start = floor(3/SI(1)) ; % 3seconds in points
for a= 1:size(GExc_hpf_Res,1) ; % changed from GExc_hpf to GExc_hpf_Res 8/4/10
    ccAltPre(a,:) = xcov(GAlt_ExcInt_hpf(a,start:prePnts),GAlt_InhInt_hpf(a,start:prePnts),'coef') ; %cross corr of alternating g exc and inh prepoints 
    ccPre(a,:) = xcov(GExc_hpf(a,start:prePnts),GInh_hpf(a,start:prePnts),'coef') ; % cross corr of single hold g exc and inh prepoints
        
    ccAltRes(a,:) = xcov(GAlt_ExcInt_hpf_Res(a,prePnts:end-postPnts),GAlt_InhInt_hpf_Res(a,prePnts:end-postPnts),'coef') ; % cross corr of alternating g exc and inh stimulus residuals
    ccRes(a,:) = xcov(GExc_hpf_Res(a,prePnts:end-postPnts),GInh_hpf_Res(a,prePnts:end-postPnts),'coef') ; % cross corr of single hold g exc and inh stimulus residuals    

    acExcPre(a,:) = xcov(GExc_hpf(a,start:prePnts),'coef') ; % auttocorrelation of prestim exc single hold
    acInhPre(a,:) = xcov(GInh_hpf(a,start:prePnts),'coef') ; % auttocorrelation of prestim inh single hold
    
    acAltExcRes(a,:) = xcov(GAlt_ExcInt_hpf_Res(a,prePnts:end-postPnts),'coef') ; % auttocorrelation of residual exc alt hold
    acAltInhRes(a,:) = xcov(GAlt_InhInt_hpf_Res(a,prePnts:end-postPnts),'coef') ; % auttocorrelation of residual inh alt hold
    
    acExcRes(a,:) = xcov(GExc_hpf_Res(a,prePnts:end-postPnts)) ; % auttocorrelation of residual exc single hold (added 8/4/10)
    acInhRes(a,:) = xcov(GInh_hpf_Res(a,prePnts:end-postPnts)) ;
    
    % cc of shuffled simultaneous recordings
    if a~=size(GExc_hpf_Res,1) ; % if your not on the last trace (% changed from GExc_hpf to GExc_hpf_Res 8/4/10)
        ccAltPre_shuf(a,:) = xcov(GAlt_ExcInt_hpf(a,start:prePnts),GAlt_InhInt_hpf(a+1,start:prePnts),'coef') ;
        ccAltRes_shuf(a,:) = xcov(GAlt_ExcInt_hpf_Res(a,prePnts:end-postPnts),GAlt_InhInt_hpf_Res(a+1,prePnts:end-postPnts),'coef') ;
    else
        ccAltPre_shuf(a,:) = xcov(GAlt_ExcInt_hpf(a,start:prePnts),GAlt_InhInt_hpf(1,start:prePnts),'coef') ;
        ccAltRes_shuf(a,:) = xcov(GAlt_ExcInt_hpf_Res(a,prePnts:end-postPnts),GAlt_InhInt_hpf_Res(1,prePnts:end-postPnts),'coef') ;
    end

end

ccAlt_Mean = xcov(GAlt_ExcInt_hpf_Mean(1,prePnts:end-postPnts),GAlt_InhInt_hpf_Mean(1,prePnts:end-postPnts),'coef') ; % cross corr of alternating g exc and inh stimulus mean
cc_Mean = xcov(GExc_hpf_Mean(1,prePnts:end-postPnts),GInh_hpf_Mean(1,prePnts:end-postPnts),'coef') ; % cross corr of single hold g exc and inh stimulus mean

% mean cross corr 
ccAltPre_Mean = mean(ccAltPre,1) ;
ccPre_Mean = mean(ccPre,1) ;
ccAltPre_shuf_Mean = mean(ccAltPre_shuf,1) ;

ccAltRes_Mean = mean(ccAltRes,1) ;
ccRes_Mean = mean(ccRes,1) ;
ccAltRes_shuf_Mean = mean(ccAltRes_shuf,1) ;

acExcPre_Mean = mean(acExcPre,1) ;
acInhPre_Mean = mean(acInhPre,1) ;

acExcRes_Mean = mean(acExcRes,1) ;
acInhRes_Mean = mean(acInhRes,1) ;

acAltExcRes_Mean = mean(acAltExcRes,1) ;
acAltInhRes_Mean = mean(acAltInhRes,1) ;

%time vectors for cc
time_ccPre = [SI(1)*([1:length(ccPre_Mean)] - (length(ccPre_Mean)+1)/2)] ; 
time_cc = [SI(1)*([1:length(cc_Mean)] - (length(cc_Mean)+1)/2)] ;

% % significance bars of cross correlation during stimulus
% numerator = max(xcov(var(GExc_hpf_Res(:,prePnts:end-postPnts),0,1),var(GInh_hpf_Res(:,prePnts:end-postPnts),0,1))...
%     +xcov(mean(GExc_hpf_Res(:,prePnts:end-postPnts),1).^2,var(GInh_hpf_Res(:,prePnts:end-postPnts),0,1))...
%     +xcov(var(GExc_hpf_Res(:,prePnts:end-postPnts),0,1),mean(GInh_hpf_Res(:,prePnts:end-postPnts),1).^2)...
%     /size(GExc_hpf_Res(:,prePnts:end-postPnts),1)) ; % taken from Brody 1999 eq 2.4
% denominator = max(acExcRes_Mean.*acInhRes_Mean) ; % corr coef normalization
% ccRes_TwoStdError = 2*sqrt(numerator/denominator) ; % 2 std of cc (if distribution is gaussian these are 95% confidence intervals)
% 
% numerator = sqrt(xcov(var(GAlt_ExcInt_hpf_Res,0,1),var(GAlt_InhInt_hpf_Res,0,1))+xcov(mean(GAlt_ExcInt_hpf_Res,1),var(GAlt_InhInt_hpf_Res,0,1))+xcov(var(GAlt_ExcInt_hpf_Res,0,1),mean(GAlt_InhInt_hpf_Res,1))/size(GExc_hpf_Res,1)) ; % taken from Brody 1999 eq 2.4
% denominator = sqrt(acAltExcRes_Mean.*acAltInhRes_Mean) ; % corr coef normalization
% ccAltRes_TwoStdError = 2*numerator./denominator ; % 2 std of cc (if distribution is gaussian these are 95% confidence intervals)


%peak and lag of cc within -200ms:+200ms
[ccAltPre_Peak,ccAltPre_Peaki] = max(abs(ccAltPre_Mean(time_ccPre<.2 & time_ccPre>-.2))) ;
ccAltPre_Peak = ccAltPre_Mean(ccAltPre_Peaki+find(time_ccPre>-.2,1)) ;

[ccPre_Peak,ccPre_Peaki] = max(abs(ccPre_Mean(time_ccPre<.2 & time_ccPre>-.2))) ;
ccPre_Peak = ccPre_Mean(ccPre_Peaki+find(time_ccPre>-.2,1)) ;

[ccAltPre_shuf_Peak,ccAltPre_shuf_Peaki] = max(abs(ccAltPre_shuf_Mean(time_ccPre<.2 & time_ccPre>-.2))) ;
ccAltPre_shuf_Peak = ccAltPre_shuf_Mean(ccAltPre_shuf_Peaki+find(time_ccPre>-.2,1)) ;

[ccAltRes_Peak,ccAltRes_Peaki] = max(abs(ccAltRes_Mean(time_cc<.2 & time_cc>-.2))) ;
ccAltRes_Peak = ccAltRes_Mean(ccAltRes_Peaki+find(time_cc>-.2,1)) ;

[ccRes_Peak,ccRes_Peaki] = max(abs(ccRes_Mean(time_cc<.2 & time_cc>-.2))) ;
ccRes_Peak = ccRes_Mean(ccRes_Peaki+find(time_cc>-.2,1)) ;

[ccAltRes_shuf_Peak,ccAltRes_shuf_Peaki] = max(abs(ccAltRes_shuf_Mean(time_cc<.2 & time_cc>-.2))) ;
ccAltRes_shuf_Peak = ccAltRes_shuf_Mean(ccAltRes_shuf_Peaki+find(time_cc>-.2,1)) ;

[ccAlt_Peak,ccAlt_Peaki] = max(ccAlt_Mean) ;
[cc_Peak,cc_Peaki] = max(cc_Mean) ;

ccAltPre_Lag = time_ccPre(ccAltPre_Peaki+find(time_ccPre>-.2,1)) ; 
ccPre_Lag= time_ccPre(ccPre_Peaki+find(time_ccPre>-.2,1)) ;
ccAltPre_shuf_Lag= time_ccPre(ccAltPre_shuf_Peaki+find(time_ccPre>-.2,1)) ;

ccAltRes_Lag= time_cc(ccAltRes_Peaki+find(time_cc>-.2,1)) ;
ccRes_Lag= time_cc(ccRes_Peaki+find(time_cc>-.2,1)) ;
ccAlt_Lag= time_cc(ccAlt_Peaki) ;
cc_Lag= time_cc(cc_Peaki) ;
ccAltRes_shuf_Lag= time_cc(ccAltRes_shuf_Peaki+find(time_cc>-.2,1)) ;


% error (mean squared) between single hold and alternating mean 
GExcError_mean = sum((GExc_hpf_Mean(prePnts:end-postPnts) - GAlt_ExcInt_hpf_Mean(prePnts:end-postPnts)).^2) ; 
GInhError_mean = sum((GInh_hpf_Mean(prePnts:end-postPnts) - GAlt_InhInt_hpf_Mean(prePnts:end-postPnts)).^2) ;

for a= 1:size(GExc_hpf,1) ; % for each trace
    GExcError(a) = sum((GExc_hpf_Mean(prePnts:end-postPnts) - GAlt_ExcInt_hpf(a,prePnts:end-postPnts)).^2) ; % difference between single hold mean and an individual alternating hold
    GInhError(a) = sum((GInh_hpf_Mean(prePnts:end-postPnts) - GAlt_InhInt_hpf_Mean(prePnts:end-postPnts)).^2) ; 

    GExcError_sh(a) = sum((GExc_hpf_Mean(prePnts:end-postPnts) - GExc_hpf(a,prePnts:end-postPnts)).^2) ; % difference between single hold mean and an individual singl hold 
    GInhError_sh(a) = sum((GInh_hpf_Mean(prePnts:end-postPnts) - GInh_hpf(a,prePnts:end-postPnts)).^2) ; 
end

% sample single hold as if it were sampled like alt hold 
GExc_rs = GExc_hpf(:,[FirstAltPnt:cyclepnts:end]) ; % resamples
GInh_rs = GInh_hpf(:,[FirstAltPnt+cyclepnts/2:cyclepnts:end]) ;

for a= 1:size(GExc_hpf,1) ; % for each trace
    GExc_rsi(a,:) = interp1([FirstAltPnt:cyclepnts:length(GExc_hpf)],GExc_rs(a,:),[1:length(GExc_hpf)],'linear','extrap') ; % resampled and interpolated
    GInh_rsi(a,:) = interp1([FirstAltPnt+cyclepnts/2:cyclepnts:length(GInh_hpf)],GInh_rs(a,:),[1:length(GInh_hpf)],'linear','extrap') ; % 
end

GExc_rsi_Mean = mean(GExc_rsi) ;
GInh_rsi_Mean = mean(GInh_rsi) ;

% get time depedant standard deviation
GAlt_ExcInt_hpf_std = std(GAlt_ExcInt_hpf) ;
GAlt_InhInt_hpf_std = std(GAlt_InhInt_hpf) ;

GExc_hpf_std = std(GExc_hpf) ;
GInh_hpf_std = std(GInh_hpf) ;

GExc_rsi_std = std(GExc_rsi) ;
GInh_rsi_std = std(GInh_rsi) ;


% variance (amplitude) of each G trace
for a= 1:size(GExc_hpf,1) ; % for each trace
    GExc_hpf_rms(a) = var(GExc_hpf(a,prePnts:end-postPnts)) ; 
    GInh_hpf_rms(a) = var(GInh_hpf(a,prePnts:end-postPnts)) ;
    
    GExc_rsi_rms(a) = var(GExc_rsi(a,prePnts:end-postPnts)) ; 
    GInh_rsi_rms(a) = var(GInh_rsi(a,prePnts:end-postPnts)) ;
    
    GAlt_ExcInt_hpf_rms(a) = var(GAlt_ExcInt_hpf(a,prePnts:end-postPnts)) ;
    GAlt_InhInt_hpf_rms(a) = var(GAlt_InhInt_hpf(a,prePnts:end-postPnts)) ;
    
    GExcPre_hpf_rms(a) = var(GExc_hpf(a,start:prePnts)) ; 
    GInhPre_hpf_rms(a) = var(GInh_hpf(a,start:prePnts)) ;
    
    GExcPre_rsi_rms(a) = var(GExc_rsi(a,start:prePnts)) ; 
    GInhPre_rsi_rms(a) = var(GInh_rsi(a,start:prePnts)) ;
    
    GAltPre_ExcInt_hpf_rms(a) = var(GAlt_ExcInt_hpf(a,start:prePnts)) ;
    GAltPre_InhInt_hpf_rms(a) = var(GAlt_InhInt_hpf(a,start:prePnts)) ;
        
end

% compare time dependant mean and std of single hold and alternating hold by getting
% the correlation coefficient and fiting a straight line to plots

[MeanExcCorr,MeanExcCorrP] = corrcoef([GExc_hpf_Mean(prePnts:end-postPnts)',GAlt_ExcInt_hpf_Mean(prePnts:end-postPnts)',GExc_rsi_Mean(prePnts:end-postPnts)']) ; % corelation coef exc sh mean and exc alt mean (3x3 matrix)
[MeanInhCorr,MeanInhCorrP] = corrcoef([GInh_hpf_Mean(prePnts:end-postPnts)',GAlt_InhInt_hpf_Mean(prePnts:end-postPnts)',GInh_rsi_Mean(prePnts:end-postPnts)']) ; % inh

[StdExcCorr,StdExcCorrP] = corrcoef([GExc_hpf_std(prePnts:end-postPnts)',GAlt_ExcInt_hpf_std(prePnts:end-postPnts)',GExc_rsi_std(prePnts:end-postPnts)']) ; % correlation coef of exc sh var and exc alt var
[StdInhCorr,StdInhCorrP] = corrcoef([GInh_hpf_std(prePnts:end-postPnts)',GAlt_InhInt_hpf_std(prePnts:end-postPnts)',GInh_rsi_std(prePnts:end-postPnts)']) ; % inh

MeanExcFitParams = nlinfit(GExc_hpf_Mean(prePnts:end-postPnts),GAlt_ExcInt_hpf_Mean(prePnts:end-postPnts),@StraightLineFit,[1,0]) ; % straight line fit comparing sh mean and alt mean
MeanInhFitParams = nlinfit(GInh_hpf_Mean(prePnts:end-postPnts),GAlt_InhInt_hpf_Mean(prePnts:end-postPnts),@StraightLineFit,[1,0]) ;

MeanExcRsiFitParams = nlinfit(GExc_rsi_Mean(prePnts:end-postPnts),GAlt_ExcInt_hpf_Mean(prePnts:end-postPnts),@StraightLineFit,[1,0]) ; % straight line fit comparing sh rsi mean and alt mean
MeanInhRsiFitParams = nlinfit(GInh_rsi_Mean(prePnts:end-postPnts),GAlt_InhInt_hpf_Mean(prePnts:end-postPnts),@StraightLineFit,[1,0]) ;

StdExcFitParams = nlinfit(GExc_hpf_std(prePnts:end-postPnts),GAlt_ExcInt_hpf_std(prePnts:end-postPnts),@StraightLineFit,[1,0]) ;
StdInhFitParams = nlinfit(GInh_hpf_std(prePnts:end-postPnts),GAlt_InhInt_hpf_std(prePnts:end-postPnts),@StraightLineFit,[1,0]) ;
   
StdExcRsiFitParams = nlinfit(GExc_rsi_std(prePnts:end-postPnts),GAlt_ExcInt_hpf_std(prePnts:end-postPnts),@StraightLineFit,[1,0]) ;
StdInhRsiFitParams = nlinfit(GInh_rsi_std(prePnts:end-postPnts),GAlt_InhInt_hpf_std(prePnts:end-postPnts),@StraightLineFit,[1,0]) ;

% power spectrums
[powerspec_xvalues, mean_powerspecExc] = PowerSpectrumFinder(GExc_hpf_Mean(:,prePnts:end-postPnts),samplerate) ; % mean single hold
[powerspec_xvalues, mean_powerspecInh] = PowerSpectrumFinder(GInh_hpf_Mean(:,prePnts:end-postPnts),samplerate) ;

unresolvable = find(powerspec_xvalues>(1/(2*cyclepnts*SI(1)))) ; % idicies of powerspec which we cannot assess with alt V (theoretical)

[powerspec_xvalues, mean_powerspecAltExc] = PowerSpectrumFinder(GAlt_ExcInt_hpf_Mean(:,prePnts:end-postPnts),samplerate) ; % mean alternanting hold
[powerspec_xvalues, mean_powerspecAltInh] = PowerSpectrumFinder(GAlt_InhInt_hpf_Mean(:,prePnts:end-postPnts),samplerate) ;
mean_powerspecAltExc(unresolvable) = 0 ; % get rid of frequencies that are unresovable
mean_powerspecAltInh(unresolvable) = 0 ; 

[powerspec_xvalues, mean_powerspecExcRes] = PowerSpectrumFinder(GExc_hpf_Res(:,prePnts:end-postPnts),samplerate) ; % residual single hold
[powerspec_xvalues, mean_powerspecInhRes] = PowerSpectrumFinder(GInh_hpf_Res(:,prePnts:end-postPnts),samplerate) ;

[powerspec_xvalues, mean_powerspecAltExcRes] = PowerSpectrumFinder(GAlt_ExcInt_hpf_Res(:,prePnts:end-postPnts),samplerate) ; % residual alt hold
[powerspec_xvalues, mean_powerspecAltInhRes] = PowerSpectrumFinder(GAlt_InhInt_hpf_Res(:,prePnts:end-postPnts),samplerate) ;
mean_powerspecAltExcRes(unresolvable) = 0 ; % get rid of frequencies that are unresovable
mean_powerspecAltInhRes(unresolvable) = 0 ;

Frac_powerspecShAltExc = mean_powerspecAltExc./mean_powerspecExc ;   % fraction of variance mean alt captures at each frequency
Frac_powerspecShAltInh = mean_powerspecAltInh./mean_powerspecInh ;

Frac_powerspecShAltExcRes = mean_powerspecAltExcRes./mean_powerspecExcRes ;   % fraction of variance residual alt captures at each frequency
Frac_powerspecShAltInhRes = mean_powerspecAltInhRes./mean_powerspecInhRes ;

FracAssessed_Exc = 1-(sum(mean_powerspecExc - mean_powerspecAltExc)/sum(mean_powerspecExc)) ; % total fraction of variance mean alt captures
FracAssessed_Inh = 1-(sum(mean_powerspecInh - mean_powerspecAltInh)/sum(mean_powerspecInh)) ;

FracAssessed_ExcRes = 1-(sum(mean_powerspecExcRes - mean_powerspecAltExcRes)/sum(mean_powerspecExcRes)) ; % total fraction of variance residual alt captures
FracAssessed_InhRes = 1-(sum(mean_powerspecInhRes - mean_powerspecAltInhRes)/sum(mean_powerspecInhRes)) ;

[powerspec_xvaluesCC, mean_powerspecAltResCC] = PowerSpectrumFinder(ccAltRes,samplerate) ; % power spectrum of alt-v cc res (correlated noise)
[powerspec_xvaluesCC, mean_powerspecAltMeanCC] = PowerSpectrumFinder(ccAlt_Mean,samplerate) ; %power spectrum of alt-v cc of mean (correlated signal)


% transfer function from single hold to alt
[LinearFilter] = LinFilterFinder(GExc_hpf(:,prePnts:end-postPnts),GAlt_ExcInt_hpf(:,prePnts:end-postPnts), samplerate, samplerate) ;
[tf_x,Exc_TransFnc] = PowerSpectrumFinder(LinearFilter,samplerate) ;

[LinearFilter] = LinFilterFinder(GInh_hpf(:,prePnts:end-postPnts),GAlt_InhInt_hpf(:,prePnts:end-postPnts), samplerate, samplerate) ;
[tf_x,Inh_TransFnc] = PowerSpectrumFinder(LinearFilter,samplerate) ;

% blanked capacitive transient in alt current
dataAltVBlanked = nans(size(dataAltV)) ;
for a = 22:51 ;        % for each data point per cycle we want to plot
    dataAltVBlanked(:,[a:cyclepnts/2:end]) = dataAltV(:,[a:cyclepnts/2:end]) ; 
end 

% figures

% figure % example alternating I exctraction
% plot(time,dataAltV(3,:))
% hold on
% plot(time,dataAltVBlanked(3,:),'k')
% plot(time,Alt_Exc(3,:),'go')
% plot(time,Alt_Inh(3,:),'ro')
% plot(time,Alt_ExcInt(3,:),'g-')
% plot(time,Alt_InhInt(3,:),'r-')
% xlabel('time (seconds)')
% ylabel('current (pA)')
% title('alternating extraction')
% 
% figure % individual example
% subplot(2,1,1)
% plot(time,GExc_hpf(3,:),'g')
% hold on
% plot(time,GAlt_ExcInt_hpf(3,:),'g--')
% title('individual comparisons')
% 
% subplot(2,1,2)
% plot(time,GInh_hpf(3,:),'r')
% hold on
% plot(time,GAlt_InhInt_hpf(3,:),'r--')
% xlabel('time (seconds)')
% ylabel('conductance (nS)')
% 
% 
% figure % individual examples of alternating hold responses
% for a=1:size(GExc_hpf,1) ;
%     subplot(size(GExc_hpf,1),1,a)
%     plot(time,GAlt_ExcInt_hpf(a,:),'g')
%     hold on
%     plot(time,GAlt_InhInt_hpf(a,:),'r')
% end
% xlabel('time (seconds)')
% ylabel('conductance (nS)')
% title('individual alternating')
% 
% figure % individual examples of single hold response
% for a=1:size(GExc_hpf,1) ;
%     subplot(size(GExc_hpf,1),1,a)
%     plot(time,GExc_hpf(a,:),'g')
%     hold on
%     plot(time,GInh_hpf(a,:),'r')
% end
% xlabel('time (seconds)')
% ylabel('conductance (nS)')
% title('individual single hold')
% 
% figure % individual examples of alternating hold residuals
% for a=1:size(GExc_hpf,1) ;
%     subplot(size(GExc_hpf,1),1,a)
%     plot(time,GAlt_ExcInt_hpf_Res(a,:),'g')
%     hold on
%     plot(time,GAlt_InhInt_hpf_Res(a,:),'r')
% end
% xlabel('time (seconds)')
% ylabel('conductance (nS)')
% title('individual alternating residuals')
% 
% figure % individual examples of single hold residuals
% for a=1:size(GExc_hpf,1) ;
%     subplot(size(GExc_hpf,1),1,a)
%     plot(time,GExc_hpf_Res(a,:),'g')
%     hold on
%     plot(time,GInh_hpf_Res(a,:),'r')
% end
% xlabel('time (seconds)')
% ylabel('conductance (nS)')
% title('individual single hold residuals')
% 
% figure % std
% subplot(1,5,1:3)
% plot(time,GExc_hpf_std,'g')
% hold on
% plot(time,GInh_hpf_std,'r')
% plot(time,GAlt_ExcInt_hpf_std,'g--')
% plot(time,GAlt_InhInt_hpf_std,'r--')
% xlabel('time (seconds)')
% ylabel('conductance (nS)')
% title('std conductances')
% 
% subplot(1,5,4)
% plot(GExc_hpf_std(prePnts:end-postPnts),GAlt_ExcInt_hpf_std(prePnts:end-postPnts),'b.')
% hold on
% plot(GExc_rsi_std(prePnts:end-postPnts),GAlt_ExcInt_hpf_std(prePnts:end-postPnts),'y.')
% plot(GExc_hpf_std(prePnts:end-postPnts),StdExcFitParams(1)*GExc_hpf_std(prePnts:end-postPnts)+StdExcFitParams(2),'b')
% plot(GExc_rsi_std(prePnts:end-postPnts),StdExcRsiFitParams(1)*GExc_rsi_std(prePnts:end-postPnts)+StdExcRsiFitParams(2),'y')
% plot([min(GExc_hpf_std(prePnts:end-postPnts)),max(GExc_hpf_std(prePnts:end-postPnts))],[min(GExc_hpf_std(prePnts:end-postPnts)),max(GExc_hpf_std(prePnts:end-postPnts))],'k-')
% text(.1,.8,['exc sh= ',num2str(StdExcCorr(1,2))],'Units','norm')
% text(.1,.7,['exc rsi= ',num2str(StdExcCorr(2,3))],'Units','norm')
% xlabel('exc sh or exc sh rsi (nS)')
% ylabel('exc alt (nS)')
% title('exc comparison')
% 
% subplot(1,5,5)
% plot(GInh_hpf_std(prePnts:end-postPnts),GAlt_InhInt_hpf_std(prePnts:end-postPnts),'b.')
% hold on
% plot(GInh_rsi_std(prePnts:end-postPnts),GAlt_InhInt_hpf_std(prePnts:end-postPnts),'y.')
% plot(GInh_hpf_std(prePnts:end-postPnts),StdInhFitParams(1)*GInh_hpf_std(prePnts:end-postPnts)+StdInhFitParams(2),'b')
% plot(GInh_rsi_std(prePnts:end-postPnts),StdInhRsiFitParams(1)*GInh_rsi_std(prePnts:end-postPnts)+StdInhRsiFitParams(2),'y')
% plot([min(GInh_hpf_std(prePnts:end-postPnts)),max(GInh_hpf_std(prePnts:end-postPnts))],[min(GInh_hpf_std(prePnts:end-postPnts)),max(GInh_hpf_std(prePnts:end-postPnts))],'k')
% text(.1,.8,['inh sh= ',num2str(StdInhCorr(1,2))],'Units','norm')
% text(.1,.7,['inh rsi= ',num2str(StdInhCorr(2,3))],'Units','norm')
% xlabel('inh sh or inh sh rsi (nS)')
% ylabel('inh alt (nS)')
% title('inh comparison')


% figure % mean traces
% subplot(1,5,1:3)
% plot(time,GExc_hpf_Mean,'g')
% hold on
% plot(time,GInh_hpf_Mean,'r')
% plot(time,GAlt_ExcInt_hpf_Mean,'g--')
% plot(time,GAlt_InhInt_hpf_Mean,'r--')
% text(3,10,['exc= ',num2str(GExcError_mean)])
% text(3,8,['inh= ',num2str(GInhError_mean)])
% xlabel('time (seconds)')
% ylabel('conductance (nS)')
% title('mean conductances')
% 
% subplot(1,5,4)
% plot(GExc_hpf_Mean(prePnts:end-postPnts),GAlt_ExcInt_hpf_Mean(prePnts:end-postPnts),'b.')
% hold on
% plot(GExc_rsi_Mean(prePnts:end-postPnts),GAlt_ExcInt_hpf_Mean(prePnts:end-postPnts),'y.')
% plot(GExc_hpf_Mean(prePnts:end-postPnts),MeanExcFitParams(1)*GExc_hpf_Mean(prePnts:end-postPnts)+MeanExcFitParams(2),'b')
% plot(GExc_rsi_Mean(prePnts:end-postPnts),MeanExcRsiFitParams(1)*GExc_rsi_Mean(prePnts:end-postPnts)+MeanExcRsiFitParams(2),'y')
% plot([min(GExc_hpf_Mean(prePnts:end-postPnts)),max(GExc_hpf_Mean(prePnts:end-postPnts))],[min(GExc_hpf_Mean(prePnts:end-postPnts)),max(GExc_hpf_Mean(prePnts:end-postPnts))],'k')
% text(.1,.8,['exc sh= ',num2str(MeanExcCorr(1,2))],'Units','norm')
% text(.1,.7,['exc rsi= ',num2str(MeanExcCorr(2,3))],'Units','norm')
% xlabel('exc sh or exc sh rsi (nS)')
% ylabel('exc alt (nS)')
% title('exc comparison')
% 
% subplot(1,5,5)
% plot(GInh_hpf_Mean(prePnts:end-postPnts),GAlt_InhInt_hpf_Mean(prePnts:end-postPnts),'b.')
% hold on
% plot(GInh_rsi_Mean(prePnts:end-postPnts),GAlt_InhInt_hpf_Mean(prePnts:end-postPnts),'y.')
% plot(GInh_hpf_Mean(prePnts:end-postPnts),MeanInhFitParams(1)*GInh_hpf_Mean(prePnts:end-postPnts)+MeanInhFitParams(2),'b')
% plot(GInh_rsi_Mean(prePnts:end-postPnts),MeanInhRsiFitParams(1)*GInh_rsi_Mean(prePnts:end-postPnts)+MeanInhRsiFitParams(2),'y')
% plot([min(GInh_hpf_Mean(prePnts:end-postPnts)),max(GInh_hpf_Mean(prePnts:end-postPnts))],[min(GInh_hpf_Mean(prePnts:end-postPnts)),max(GInh_hpf_Mean(prePnts:end-postPnts))],'k')
% text(.1,.8,['inh sh= ',num2str(MeanInhCorr(1,2))],'Units','norm')
% text(.1,.7,['inh rsi= ',num2str(MeanInhCorr(2,3))],'Units','norm')
% xlabel('inh sh or inh sh rsi (nS)')
% ylabel('inh alt (nS)')
% title('inh comparison')

% figure % mean traces using resampled single hold
% plot(time,GExc_rsi_Mean,'g')
% hold on
% plot(time,GInh_rsi_Mean,'r')
% plot(time,GAlt_ExcInt_hpf_Mean,'g--')
% plot(time,GAlt_InhInt_hpf_Mean,'r--')
% text(3,10,['exc= ',num2str(GExcError_mean)])
% text(3,8,['inh= ',num2str(GInhError_mean)])
% xlabel('time (seconds)')
% ylabel('conductance (nS)')
% title('mean conductances (resampled)')
% 
% 
% figure % rms (amplitude) of Gtraces
% plot(1,GExc_hpf_rms,'g*')
% hold on
% plot(1,GInh_hpf_rms,'r*')
% plot(2,GExc_rsi_rms,'g*')
% plot(2,GInh_rsi_rms,'r*')
% plot(3,GAlt_ExcInt_hpf_rms,'g*')
% plot(3,GAlt_InhInt_hpf_rms,'r*')
% a = gca ;
% set(a,'XTick',[1:3]) ;
% set(a,'XTickLabel',['sh    ';'sh rsi';'alt   ']) ;
% 
% figure % cross correlations
% plot(time_ccPre,ccAltPre_Mean)
% hold on
% plot(time_ccPre,ccPre_Mean,'r')
% plot(time_cc,ccAltRes_Mean,'g')
% plot(time_cc,ccRes_Mean,'y')
% plot(time_cc,ccAlt_Mean,'k')
% plot(time_cc,cc_Mean,'c')
% plot(time_ccPre,ccAltPre_shuf_Mean,'b--')
% plot(time_cc,ccAltRes_shuf_Mean,'g--')
% plot(time_ccPre,acExcPre_Mean,'b:')
% plot(time_ccPre,acInhPre_Mean,'r:')
% plot(time_cc,acAltExcRes_Mean,'b.-')
% plot(time_cc,acAltInhRes_Mean,'r.-')
% legend('alternating prestim','single hold prestim','alternating residual','single hold residuals','alternating meanG stimulus','single hold meanG stimulus','alt prestim shuffled','alt resid shuffled','autocorr sh exc','autocorr sh inh','acAltExcRes_Mean','acAltInhRes_Mean') 
% title('cross correlations')
% xlabel('time (seconds)')
% ylabel('fraction shared')
% 
% plot(ccAltPre_Lag,ccAltPre_Peak,'b*') 
% plot(ccPre_Lag,ccPre_Peak,'r*') 
% plot(ccAltRes_Lag,ccAltRes_Peak,'g*') 
% plot(ccRes_Lag,ccRes_Peak,'y*')
% plot(ccAlt_Lag,ccAlt_Peak,'k*')
% plot(cc_Lag,cc_Peak,'c*')
% 
% figure % individual cross correlations
% subplot(2,2,1)
% plot(time_ccPre,ccAltPre,'b')
% hold on
% plot(time_ccPre, ccAltPre_Mean,'b','Linewidth',2)
% plot(time_ccPre,ccPre,'r')
% plot(time_ccPre, ccPre_Mean,'r','Linewidth',2)
% axis([-.1 .1 -.5 1])
% title('preStim')
% xlabel('time (seconds)')
% ylabel('fraction shared')
% 
% subplot(2,2,2)
% plot(time_cc,ccAltRes,'b')
% hold on
% plot(time_cc, ccAltRes_Mean,'b','Linewidth',2)
% plot(time_cc,ccRes,'r')
% plot(time_cc, ccRes_Mean,'r','Linewidth',2)
% axis([-.1 .1 -.5 1])
% title('Residuals')
% xlabel('time (seconds)')
% ylabel('fraction shared')
% 
% subplot(2,2,3)
% plot(time_cc, ccAlt_Mean,'b','Linewidth',2)
% hold on
% plot(time_cc, cc_Mean,'r','Linewidth',2)
% axis([-.1 .1 -.5 1])
% title('Stimulus mean')
% xlabel('time (seconds)')
% ylabel('fraction shared')
% 
% subplot(2,2,4)
% plot(time_ccPre,ccAltPre_shuf,'y')
% hold on
% plot(time_ccPre,ccAltPre_shuf_Mean,'y','Linewidth',2)
% plot(time_cc,ccAltRes_shuf,'c')
% plot(time_cc, ccAltRes_shuf_Mean,'c','Linewidth',2)
% axis([-.1 .1 -.5 1])
% title('Shuffled simultaneous')
% legend('prestim','','Residual','')
% xlabel('time (seconds)')
% ylabel('fraction shared')

% figure % power spectrum of mean traces
% subplot(2,2,1)
% plot(powerspec_xvalues, mean_powerspecExc,'g')
% hold on
% plot(powerspec_xvalues, mean_powerspecInh,'r')
% plot(powerspec_xvalues, mean_powerspecAltExc,'b')
% plot(powerspec_xvalues, mean_powerspecAltInh,'y')
% h = gca ;
% set(h,'XScale','log','YScale','log')
% legend('exc sh', 'inh sh', 'exc alt', 'inh alt')
% title('mean signal')
% 
% subplot(2,2,2) % transfer functions 
% plot(tf_x,Exc_TransFnc,'g')
% hold on
% plot(tf_x,Inh_TransFnc,'r')
% h = gca ;
% set(h,'XScale','log','YScale','log')
% xlabel('frequency')
% ylabel('gain')
% title('transfer functions')
% 
% subplot(2,2,3) % power spectrum of residuals
% plot(powerspec_xvalues, mean_powerspecExcRes,'g')
% hold on
% plot(powerspec_xvalues, mean_powerspecInhRes,'r')
% plot(powerspec_xvalues, mean_powerspecAltExcRes,'b')
% plot(powerspec_xvalues, mean_powerspecAltInhRes,'y')
% h = gca ;
% set(h,'XScale','log','YScale','log')
% legend('exc sh', 'inh sh', 'exc alt', 'inh alt')
% title('residuals')
% 
% subplot(2,2,4) % power spectrum of residuals cc and mean cc from alt
% plot(powerspec_xvaluesCC, mean_powerspecAltResCC,'b')
% hold on
% plot(powerspec_xvaluesCC, mean_powerspecAltMeanCC,'r')
% h = gca ;
% set(h,'XScale','log','YScale','log')
% legend('cc residuals alt', 'cc mean alt')
% title('alt cc spectrum')
% 
% subplot(2,2,4) % difference of power spectra
% plot(powerspec_xvalues,Frac_powerspecShAltExcRes,'g')
% hold on
% plot(powerspec_xvalues,Frac_powerspecShAltInhRes,'r')
% h = gca ;
% set(h,'XScale','log','YScale','linear')
% legend(['exc',num2str(FracAssessed_ExcRes)],['inh',num2str(FracAssessed_InhRes) ])

% figure
% plot(GExc_hpf_dot,'g--')
% hold on
% plot(GInh_hpf_dot,'r--')
% plot(GAlt_ExcInt_hpf_dot,'g')
% plot(GAlt_InhInt_hpf_dot,'r')

%FOR IGOR

% ForIgor.nada = [] ;

% % comparing mean single hold vs mean alt g
% identifier = ['LineCorrCoefExc',id,num2str(A)] ;
% ForIgor.(identifier) = MeanExcCorr(1,2) ; 
% 
% identifier = ['LineCorrCoefInh',id,num2str(A)] ;
% ForIgor.(identifier) = MeanInhCorr(1,2) ; 
% 
% identifier = ['LineSlopeExc',id,num2str(A)] ;
% ForIgor.(identifier) = MeanExcFitParams(1) ;
% 
% identifier = ['LineSlopeInh',id,num2str(A)] ;
% ForIgor.(identifier) = MeanInhFitParams(1) ;
% 
% % comparing time dependant std single hold vs mean alt g
% identifier = ['LineStdCorrCoefExc',id,num2str(A)] ;
% ForIgor.(identifier) = StdExcCorr(1,2) ; 
% 
% identifier = ['LineStdCorrCoefInh',id,num2str(A)] ;
% ForIgor.(identifier) = StdInhCorr(1,2) ; 
% 
% identifier = ['LineStdSlopeExc',id,num2str(A)] ;
% ForIgor.(identifier) = StdExcFitParams(1) ; 
% 
% identifier = ['LineStdSlopeInh',id,num2str(A)] ;
% ForIgor.(identifier) = StdInhFitParams(1) ;
% 
% % comparing time dependant std resampled single hold vs alt g
% identifier = ['LineStdCorrCoefExcRsi',id,num2str(A)] ;
% ForIgor.(identifier) = StdExcCorr(2,3) ; 
% 
% identifier = ['LineStdCorrCoefInhRsi',id,num2str(A)] ;
% ForIgor.(identifier) = StdInhCorr(2,3) ; 
% 
% identifier = ['LineStdSlopeExcRsi',id,num2str(A)] ;
% ForIgor.(identifier) = StdExcRsiFitParams(1) ; 
% 
% identifier = ['LineStdSlopeInhRsi',id,num2str(A)] ;
% ForIgor.(identifier) = StdInhRsiFitParams(1) ;
% 
% % variance captured
% identifier = ['FracCapExc',id,num2str(A)] ;
% ForIgor.(identifier) = FracAssessed_Exc ;
% 
% identifier = ['FracCapInh',id,num2str(A)] ;
% ForIgor.(identifier) = FracAssessed_Inh ;
% 
% identifier = ['FracCapExcRes',id,num2str(A)] ;
% ForIgor.(identifier) = FracAssessed_ExcRes ;
% 
% identifier = ['FracCapInhRes',id,num2str(A)] ;
% ForIgor.(identifier) = FracAssessed_InhRes ;

% % mean cc
% identifier = ['ccPretime',id,num2str(A)] ;
% ForIgor.(identifier) = time_ccPre ;
% 
% identifier = ['ccAltPreMean',id,num2str(A)] ;
% ForIgor.(identifier) = ccAltPre_Mean ;
% 
% identifier = ['ccPreMean',id,num2str(A)] ;
% ForIgor.(identifier) = ccPre_Mean ;
% 
% identifier = ['ccShufPreMean',id,num2str(A)] ;
% ForIgor.(identifier) = ccAltPre_shuf_Mean ;
% 
identifier = ['ccRestime',id,num2str(A)] ;
ForIgor.(identifier) = time_cc ;

identifier = ['ccAltResMean',id,num2str(A)] ;
ForIgor.(identifier) = ccAltRes_Mean ;

identifier = ['ccResMean',id,num2str(A)] ;
ForIgor.(identifier) = ccRes_Mean ;
% 
% identifier = ['ccShufResMean',id,num2str(A)] ;
% ForIgor.(identifier) = ccAltRes_shuf_Mean ;
% 
% % mean cc lags and peaks
% identifier = ['ccAltPreLag',id,num2str(A)] ;
% ForIgor.(identifier) = ccAltPre_Lag ;
% 
% identifier = ['ccPreLag',id,num2str(A)] ;
% ForIgor.(identifier) = ccPre_Lag ;
% 
% identifier = ['ccShufPreLag',id,num2str(A)] ;
% ForIgor.(identifier) = ccAltPre_shuf_Lag ;
% 
% identifier = ['ccAltResLag',id,num2str(A)] ;
% ForIgor.(identifier) = ccAltRes_Lag ;
% 
% identifier = ['ccResLag',id,num2str(A)] ;
% ForIgor.(identifier) = ccRes_Lag ;
% 
% identifier = ['ccShufResLag',id,num2str(A)] ;
% ForIgor.(identifier) = ccAltRes_shuf_Lag ;
% 
% identifier = ['ccAltPrePeak',id,num2str(A)] ;
% ForIgor.(identifier) = ccAltPre_Peak ;
% 
% identifier = ['ccPrePeak',id,num2str(A)] ;
% ForIgor.(identifier) = ccPre_Peak ;
% 
% identifier = ['ccShufPrePeak',id,num2str(A)] ;
% ForIgor.(identifier) = ccAltPre_shuf_Peak ;
%
identifier = ['ccAltResPeak',id,num2str(A)] ;
ForIgor.(identifier) = ccAltRes_Peak ;

identifier = ['ccResPeak',id,num2str(A)] ;
ForIgor.(identifier) = ccRes_Peak ;
% 
% identifier = ['ccShufResPeak',id,num2str(A)] ;
% ForIgor.(identifier) = ccAltRes_shuf_Peak ;
% 
% % mean correlation coeficient (mean cc at time lag 0)
% identifier = ['corrCoefAltPreMean',id,num2str(A)] ;
% ForIgor.(identifier) = ccAltPre_Mean(time_ccPre==0) ;
% 
% identifier = ['corrCoefPreMean',id,num2str(A)] ;
% ForIgor.(identifier) = ccPre_Mean(time_ccPre==0) ;
% 
% identifier = ['corrCoefShufPreMean',id,num2str(A)] ;
% ForIgor.(identifier) = ccAltPre_shuf_Mean(time_ccPre==0) ;
% 
% identifier = ['corrCoefAltResMean',id,num2str(A)] ;
% ForIgor.(identifier) = ccAltRes_Mean(time_cc==0) ;
% 
% identifier = ['corrCoefResMean',id,num2str(A)] ;
% ForIgor.(identifier) = ccRes_Mean(time_cc==0) ;
% 
% identifier = ['corrCoefShufResMean',id,num2str(A)] ;
% ForIgor.(identifier) = ccAltRes_shuf_Mean(time_cc==0) ;

% individual cc
% for a=1:size(ccAltPre,1) ; 
%     identifier = ['ccAltPre',num2str(a),id,num2str(A)] ;
%     ForIgor.(identifier) = ccAltPre(a,:) ;
%     
%     identifier = ['ccPre',num2str(a),id,num2str(A)] ;
%     ForIgor.(identifier) = ccPre(a,:) ;    
%     
%     identifier = ['ccShufPre',num2str(a),id,num2str(A)] ;
%     ForIgor.(identifier) = ccAltPre_shuf(a,:) ;
%     
%     identifier = ['ccAltRes',num2str(a),id,num2str(A)] ;
%     ForIgor.(identifier) = ccAltRes(a,:) ;
%     
%     identifier = ['ccRes',num2str(a),id,num2str(A)] ;
%     ForIgor.(identifier) = ccRes(a,:) ;    
%     
%     identifier = ['ccShufRes',num2str(a),id,num2str(A)] ;
%     ForIgor.(identifier) = ccAltRes_shuf(a,:) ;
% end
% 
% % individual alt and sh conductances and residuals
% for a=1:size(GExc_hpf,1) ; 
%     
%     identifier = ['exGshExc',num2str(a),id,num2str(A)] ;
%     ForIgor.(identifier) = GExc_hpf(a,:) ;
% 
%     identifier = ['exGshInh',num2str(a),id,num2str(A)] ;
%     ForIgor.(identifier) = GInh_hpf(a,:) ;    
%     
%     identifier = ['exGaltExc',num2str(a),id,num2str(A)] ;
%     ForIgor.(identifier) = GAlt_ExcInt_hpf(a,:) ;
% 
%     identifier = ['exGaltInh',num2str(a),id,num2str(A)] ;
%     ForIgor.(identifier) = GAlt_InhInt_hpf(a,:) ;
%     
%     identifier = ['exGaltExcRes',num2str(a),id,num2str(A)] ;
%     ForIgor.(identifier) = GAlt_ExcInt_hpf_Res(a,:) ;
% 
%     identifier = ['exGaltInhRes',num2str(a),id,num2str(A)] ;
%     ForIgor.(identifier) = GAlt_InhInt_hpf_Res(a,:) ;   
%  
%     identifier = ['exGshRsExc',num2str(a),id,num2str(A)] ;
%     ForIgor.(identifier) = GExc_rsi(a,:) ;
% 
%     identifier = ['exGshRsInh',num2str(a),id,num2str(A)] ;
%     ForIgor.(identifier) = GInh_rsi(a,:) ; 
% 
% end
% 
% % mean alt and sh conductances
% 
% identifier = ['meanGshExc',id,num2str(A)] ;
% ForIgor.(identifier) = GExc_hpf_Mean ;
% 
% identifier = ['meanGshInh',id,num2str(A)] ;
% ForIgor.(identifier) = GInh_hpf_Mean ;
% 
% identifier = ['meanGaltExc',id,num2str(A)] ;
% ForIgor.(identifier) = GAlt_ExcInt_hpf_Mean ;
% 
% identifier = ['meanGaltInh',id,num2str(A)] ;
% ForIgor.(identifier) = GAlt_InhInt_hpf_Mean ;   
% 
% identifier = ['meanGshRsExc',id,num2str(A)] ;
% ForIgor.(identifier) = GExc_rsi_Mean ;
% 
% identifier = ['meanGshRsInh',id,num2str(A)] ;
% ForIgor.(identifier) = GInh_rsi_Mean ;
% 
% % transfer functions from sh to alt individs
% 
% identifier = ['TransFncExc',id,num2str(A)] ;
% ForIgor.(identifier) = Exc_TransFnc ;
% 
% identifier = ['TransFncInh',id,num2str(A)] ;
% ForIgor.(identifier) = Inh_TransFnc ;
% 
% identifier = ['TransFncX',id,num2str(A)] ;
% ForIgor.(identifier) = tf_x ;
% 
% % voltage command and light stimulus
% 
% identifier = ['VoltageCommand',id,num2str(A)] ;
% ForIgor.(identifier) = voltageCommand(3,:) ;
% 
% identifier = ['LightCommand',id,num2str(A)] ;
% ForIgor.(identifier) = lightCommand(3,:) ;
% 
% % current alt voltage trial
% 
% identifier = ['IaltV',id,num2str(A)] ;
% ForIgor.(identifier) = dataAltV(3,:) ;
% 
% identifier = ['IaltVblanked',id,num2str(A)] ;
% ForIgor.(identifier) = dataAltVBlanked(3,:) ;
% 
% identifier = ['IaltVExc',id,num2str(A)] ;
% ForIgor.(identifier) = Alt_Exc(3,:) ;
% 
% identifier = ['IaltVInh',id,num2str(A)] ;
% ForIgor.(identifier) = Alt_Inh(3,:) ;
% 
% identifier = ['IaltVExcInt',id,num2str(A)] ;
% ForIgor.(identifier) = Alt_ExcInt(3,:) ;
% 
% identifier = ['IaltVInhInt',id,num2str(A)] ;
% ForIgor.(identifier) = Alt_InhInt(3,:) ;
% 
% % time
% identifier = ['time',id,num2str(A)] ;
% ForIgor.(identifier) = time ;
% 
% identifier = ['shExcdot',id,num2str(A)] ;
% ForIgor.(identifier) = GExc_hpf_dot ;
% 
% identifier = ['shInhdot',id,num2str(A)] ;
% ForIgor.(identifier) = GInh_hpf_dot ;
% 
% identifier = ['altExcdot',id,num2str(A)] ;
% ForIgor.(identifier) = GAlt_ExcInt_hpf_dot ;
% 
% identifier = ['altInhdot',id,num2str(A)] ;
% ForIgor.(identifier) = GAlt_InhInt_hpf_dot ;

% % residual option cc data
% if residualOption~=0 ;
%     identifier = ['ccAltResPeakResOpt',id,num2str(A)] ;
%     ForIgor.(identifier) = ccAltRes_Peak ;
% 
%     identifier = ['ccResPeakResOpt',id,num2str(A)] ;
%     ForIgor.(identifier) = ccRes_Peak ;
% 
%     identifier = ['ccAltResMeanResOpt',id,num2str(A)] ;
%     ForIgor.(identifier) = ccAltRes_Mean ;
% 
%     identifier = ['ccResMeanResOpt',id,num2str(A)] ;
%     ForIgor.(identifier) = ccRes_Mean ;
% end



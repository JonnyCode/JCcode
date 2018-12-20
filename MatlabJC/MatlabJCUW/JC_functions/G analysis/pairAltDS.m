function ForIgor = pairAltDS(Input,Parameters,id,A) ; 

% this function will analyze data from pairs of ds cells recroded in alt
% voltage exp
% JC 4/26/10
% JC 1/20/11 added MonitorAnalysis code to more accurately extract framerate 

% get data
[fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/mouse/',Input(A).cellname]);

epochs = str2num(Input(A).(id)) ;
round = 0 ;
for a = 1:3:length(epochs) ; % for each spike epoch
    round = round +1 ;

    [voltageCommand1(round,:), error] = ITCReadEpochStm(epochs(a+2), 3,fp) ; % get voltage command
    [voltageCommand2(round,:), error] = ITCReadEpochStm(epochs(a+2), 4,fp) ; 
    
    [dataAltV1(round,:), error] = ITCReadEpoch(epochs(a+2), 0, fp) ;    % get data
    [dataAltV2(round,:), error] = ITCReadEpoch(epochs(a+2), 1, fp) ;
    
    
    if voltageCommand1(round,1)==min(voltageCommand1(round,1:500)) ;
        [dataExc1(round,:), error] = ITCReadEpoch(epochs(a), 0, fp) ;    %#ok<*AGROW> % get data
        [dataInh1(round,:), error] = ITCReadEpoch(epochs(a+1), 0, fp) ;    % get data
    else
        [dataExc1(round,:), error] = ITCReadEpoch(epochs(a+1), 0, fp) ;    %#ok<*AGROW> % get data
        [dataInh1(round,:), error] = ITCReadEpoch(epochs(a), 0, fp) ;    % get data
    end
        
    if voltageCommand2(round,1)==min(voltageCommand2(round,1:500)) ;
        [dataExc2(round,:), error] = ITCReadEpoch(epochs(a), 1, fp) ;    %#ok<*AGROW> % get data
        [dataInh2(round,:), error] = ITCReadEpoch(epochs(a+1), 1, fp) ;    % get data
    else
        [dataExc2(round,:), error] = ITCReadEpoch(epochs(a+1), 1, fp) ;    %#ok<*AGROW> % get data
        [dataInh2(round,:), error] = ITCReadEpoch(epochs(a), 1, fp) ;    % get data
    end      

    
    [SI(round), error] = ITCGetSamplingInterval(epochs(a+2), fp); % get sampling interval
    SI(round) = SI(round) * 1e-6; % Sampling interval in sec
end

voltageCommand1 = voltageCommand1(:,1:length(dataExc1)) ;
voltageCommand2 = voltageCommand2(:,1:length(dataExc1)) ;

if Input(A).ITC18flag == 1 ;
    SI = SI*1.25 ;
end

SpatialStimParams = hdf5load(['~/Data/mouse/',Input(A).cellname,'_spatial.h5']) ; % load spatial stim params

for a = 1:length(epochs) ;
    [Monitor(a,:), error] = ITCReadEpoch(epochs(a), 2, fp) ;    % get monitor data
end
[framerate_mean,framerate_range,startAmp] = MonitorAnalysis(Monitor, length(epochs), SI) ; % get frame rate
frameRate = mean(framerate_mean) ;
%frameRate = 60 ;

for a = 1:length(epochs) ;
    StrucString = ['params_epoch',num2str(epochs(a))] ; 
    Struct = SpatialStimParams.(StrucString) ;

    BarAngles(a,:) = Struct.BarAngle ;
    BarWidth(a) = Struct.BarWidth ;
    BarSpeed(a) = Struct.BarSpeed ;
    
    prePnts(a) = floor(Struct.spatial_prepts/(frameRate*SI(1))) ;
    postPnts(a) = floor(Struct.spatial_postpts/(frameRate*SI(1))) ;
    groupPnts(a) = floor((length(dataExc1)-prePnts(a)-postPnts(a))/length(BarAngles(a,:))) ;
    
    NumBars(a) = length(BarAngles(a,:)) ; % number of bars shown per trial
    BarPnts(a) = floor((length(dataExc1)-prePnts(a)-postPnts(a))/NumBars(a)) ; % number of points the bar is presented + interbar points
    OnPnts(a) = floor((BarWidth(a)/BarSpeed(a))/(frameRate*SI(1))) ; % number of points before end of bar appears triggering off response
end

cyclepnts = 100 ; % number of sample points in a cycle, pnts between leaving hold1 and returning (Also gets rid of first cycle)    
FirstAltPnt = (cyclepnts/2)+1 ; % first sample point you want to plot after begining of step from alternation 
LastAltPnt = FirstAltPnt ;  % last "                                                                    "

samplerate = 1/SI(1) ; % Hz at which data was collected
time = [SI(1):SI(1):SI(1)*length(dataExc1(1,:))] ; % time vector in seconds

% WINDOWING ALTERNATING V DATA
Alt_Exc1 = NaN(size(dataAltV1)) ; % make a vector of NaNs that will serve as base for ploting alternating coductances at hold1
Alt_Inh1 = NaN(size(dataAltV1)) ; % make a vector of NaNs that will serve as base for ploting alternating coductances at hold2

Alt_Exc2 = NaN(size(dataAltV2)) ; % make a vector of NaNs that will serve as base for ploting alternating coductances at hold1
Alt_Inh2 = NaN(size(dataAltV2)) ; % make a vector of NaNs that will serve as base for ploting alternating coductances at hold2

for a = FirstAltPnt:LastAltPnt ;        % for each data point per cycle we want to plot  
    if voltageCommand1(1,1)==min(voltageCommand1(1,1:500)) ; 
        Alt_Exc1(:,[a:cyclepnts:end]) = dataAltV1(:,[a:cyclepnts:end]) ; 
        Alt_Inh1(:,[a+cyclepnts/2:cyclepnts:end]) = dataAltV1(:,[a+cyclepnts/2:cyclepnts:end]) ; 
    else
        Alt_Inh1(:,[a:cyclepnts:end]) = dataAltV1(:,[a:cyclepnts:end]) ; 
        Alt_Exc1(:,[a+cyclepnts/2:cyclepnts:end]) = dataAltV1(:,[a+cyclepnts/2:cyclepnts:end]) ;  
    end

    if voltageCommand2(1,1)==min(voltageCommand2(1,1:500)) ; 
        Alt_Exc2(:,[a:cyclepnts:end]) = dataAltV2(:,[a:cyclepnts:end]) ; 
        Alt_Inh2(:,[a+cyclepnts/2:cyclepnts:end]) = dataAltV2(:,[a+cyclepnts/2:cyclepnts:end]) ; 
    else
        Alt_Inh2(:,[a:cyclepnts:end]) = dataAltV2(:,[a:cyclepnts:end]) ; 
        Alt_Exc2(:,[a+cyclepnts/2:cyclepnts:end]) = dataAltV2(:,[a+cyclepnts/2:cyclepnts:end]) ;
    end
end 

numTrials = length(epochs)/3 ;
% interpolate alternating current data
for a = 1:numTrials ; % for each trial
    b = find(isnan(Alt_Exc1(a,:)) == 0) ;     % find all the indices that are not nans
    Alt_ExcInt1(a,:) = interp1(b,Alt_Exc1(a,b),[1:length(Alt_Exc1)],'linear','extrap') ; % interpolate to find values that were not sampled
    clear b
    
    b = find(isnan(Alt_Exc2(a,:)) == 0) ;     % find all the indices that are not nans
    Alt_ExcInt2(a,:) = interp1(b,Alt_Exc2(a,b),[1:length(Alt_Exc2)],'linear','extrap') ; % interpolate to find values that were not sampled
    clear b
    
    b = find(isnan(Alt_Inh1(a,:)) == 0) ;     % find all the indices that are not nans
    Alt_InhInt1(a,:) = interp1(b,Alt_Inh1(a,b),[1:length(Alt_Inh1)],'linear','extrap') ; % interpolate to find values that were not sampled
    clear b
    
    b = find(isnan(Alt_Inh2(a,:)) == 0) ;     % find all the indices that are not nans
    Alt_InhInt2(a,:) = interp1(b,Alt_Inh2(a,b),[1:length(Alt_Inh2)],'linear','extrap') ; % interpolate to find values that were not sampled
    clear b
end

% low pass filter and remove electrical crap
dataExc_lpf1 = lowPassFilter(dataExc1, samplerate, 5000) ; %(signal,samplerate,cutoff frequ (hz))
dataInh_lpf1 = lowPassFilter(dataInh1, samplerate, 5000) ; 
Alt_ExcInt_lpf1 = lowPassFilter(Alt_ExcInt1, samplerate, 5000) ;
Alt_InhInt_lpf1 = lowPassFilter(Alt_InhInt1, samplerate, 5000) ;

dataExc_lpf2 = lowPassFilter(dataExc2, samplerate, 5000) ; %(signal,samplerate,cutoff frequ (hz))
dataInh_lpf2 = lowPassFilter(dataInh2, samplerate, 5000) ; 
Alt_ExcInt_lpf2 = lowPassFilter(Alt_ExcInt2, samplerate, 5000) ;
Alt_InhInt_lpf2 = lowPassFilter(Alt_InhInt2, samplerate, 5000) ;


% high pass filter to remove slow drift
% greg S. sent this to me to help implement a butterworth and avoid ringing
F=1 ; % filter cuttoff
Wn = F*SI(1); %normalized frequency cutoff
[z, p, k] = butter(1,Wn,'high'); %
[sos,g]=zp2sos(z,p,k); 
myfilt=dfilt.df2sos(sos,g);

dataExc_hpf1 = filter(myfilt,dataExc_lpf1')'; % filter implementation
dataInh_hpf1 = filter(myfilt,dataInh_lpf1')'; 
Alt_ExcInt_hpf1 = filter(myfilt,Alt_ExcInt_lpf1')'; 
Alt_InhInt_hpf1 = filter(myfilt,Alt_InhInt_lpf1')';

dataExc_hpf2 = filter(myfilt,dataExc_lpf2')'; % filter implementation
dataInh_hpf2 = filter(myfilt,dataInh_lpf2')'; 
Alt_ExcInt_hpf2 = filter(myfilt,Alt_ExcInt_lpf2')'; 
Alt_InhInt_hpf2 = filter(myfilt,Alt_InhInt_lpf2')';

% change to conductances and subtract off means 
for a = 1:numTrials ; % for each trial
    GExc_hpf1(a,:) = dataExc_hpf1(a,:)/-61 - mean(dataExc_hpf1(a,prePnts:end-postPnts)/-61) ; % get conductance from stable currrents
    GInh_hpf1(a,:) = dataInh_hpf1(a,:)/61 - mean(dataInh_hpf1(a,prePnts:end-postPnts)/61) ; 
    GAlt_ExcInt_hpf1(a,:) = Alt_ExcInt_hpf1(a,:)/-61 - mean(Alt_ExcInt_hpf1(a,prePnts:end-postPnts)/-61) ; % get conductance from alt current
    GAlt_InhInt_hpf1(a,:) = Alt_InhInt_hpf1(a,:)/61 - mean(Alt_InhInt_hpf1(a,prePnts:end-postPnts)/61) ; %#ok<*AGROW>
    
    GExc_hpf2(a,:) = dataExc_hpf2(a,:)/-61 - mean(dataExc_hpf2(a,prePnts:end-postPnts)/-61) ; % get conductance from stable currrents
    GInh_hpf2(a,:) = dataInh_hpf2(a,:)/61 - mean(dataInh_hpf2(a,prePnts:end-postPnts)/61) ; 
    GAlt_ExcInt_hpf2(a,:) = Alt_ExcInt_hpf2(a,:)/-61 - mean(Alt_ExcInt_hpf2(a,prePnts:end-postPnts)/-61) ; % get conductance from alt current
    GAlt_InhInt_hpf2(a,:) = Alt_InhInt_hpf2(a,:)/61 - mean(Alt_InhInt_hpf2(a,prePnts:end-postPnts)/61) ; %#ok<*AGROW>
end

offsetAltExc1 = min(min(GAlt_ExcInt_hpf1(:,prePnts:end-postPnts))) ;
offsetExc1 = min(min(GExc_hpf1(:,prePnts:end-postPnts))) ;

offsetAltInh1 = min(min(GAlt_InhInt_hpf1(:,prePnts:end-postPnts))) ;
offsetInh1 = min(min(GInh_hpf1(:,prePnts:end-postPnts))) ;

offsetAltExc2 = min(min(GAlt_ExcInt_hpf2(:,prePnts:end-postPnts))) ;
offsetExc2 = min(min(GExc_hpf2(:,prePnts:end-postPnts))) ;

offsetAltInh2 = min(min(GAlt_InhInt_hpf2(:,prePnts:end-postPnts))) ;
offsetInh2 = min(min(GInh_hpf2(:,prePnts:end-postPnts))) ;

% ofset G (assumes all g have same mean and 1 min) 
for a = 1:numTrials ; % for each trial
    GExc_hpf1(a,:) = GExc_hpf1(a,:) - offsetExc1 ; % offsets
    GInh_hpf1(a,:) = GInh_hpf1(a,:) - offsetInh1 ; 
    GAlt_ExcInt_hpf1(a,:) = GAlt_ExcInt_hpf1(a,:) - offsetAltExc1 ; 
    GAlt_InhInt_hpf1(a,:) = GAlt_InhInt_hpf1(a,:) - offsetAltInh1 ;
    
    GExc_hpf2(a,:) = GExc_hpf2(a,:) - offsetExc2 ; % offsets
    GInh_hpf2(a,:) = GInh_hpf2(a,:) - offsetInh2 ; 
    GAlt_ExcInt_hpf2(a,:) = GAlt_ExcInt_hpf2(a,:) - offsetAltExc2 ; 
    GAlt_InhInt_hpf2(a,:) = GAlt_InhInt_hpf2(a,:) - offsetAltInh2 ;
end

figure % example alternating I exctraction
subplot(2,1,1)
plot(time,dataAltV1(3,:))
hold on
plot(time,Alt_Exc1(3,:),'go')
plot(time,Alt_Inh1(3,:),'ro')
plot(time,Alt_ExcInt1(3,:),'g-')
plot(time,Alt_InhInt1(3,:),'r-')
xlabel('time (seconds)')
ylabel('current (pA)')
title('alternating extraction')

subplot(2,1,2)
plot(time,dataAltV2(3,:))
hold on
plot(time,Alt_Exc2(3,:),'go')
plot(time,Alt_Inh2(3,:),'ro')
plot(time,Alt_ExcInt2(3,:),'g-')
plot(time,Alt_InhInt2(3,:),'r-')
xlabel('time (seconds)')
ylabel('current (pA)')
title('alternating extraction')

% to save memmory clear unused variables
clearvars -except GExc_hpf1 GInh_hpf1 GAlt_ExcInt_hpf1 GAlt_InhInt_hpf1...
    GExc_hpf2 GInh_hpf2 GAlt_ExcInt_hpf2 GAlt_InhInt_hpf2...
    prePnts postPnts BarPnts OnPnts BarAngles NumBars SI Monitor A Input fp...
    frameRate SpatialStimParams id time numTrials epochs

UniqueBars = unique(BarAngles(:)) ; % unique bar angles
NumUniqueBars = length(UniqueBars) ; % number of unique bar angles

for a=1:NumUniqueBars ;
    gExc1{a} = [] ;
    gInh1{a} = [] ;
    gAlt_Exc1{a} = [] ;
    gAlt_Inh1{a} = [] ;

    gExc2{a} = [] ;
    gInh2{a} = [] ; 
    gAlt_Exc2{a} = [] ;
    gAlt_Inh2{a} = [] ;
end

% cut up and arrange data in array by bar direction and on/off responses
r=0;
[Ba,i] = sort(BarAngles,2) ; % sort each row of bar angle in accending order
for b= 1:3:size(BarAngles,1) ; % on each trial
    r=r+1 ;
    for a = 1:NumBars(b) ; % for each bar
        Ubar = find(UniqueBars==Ba(b,a)) ;
        
        gExc1{Ubar} = [gExc1{Ubar}; GExc_hpf1(r,prePnts(b)+BarPnts(b)*(i(b,a)-1):prePnts(b)+BarPnts(b)*i(b,a))] ;
        gInh1{Ubar} = [gInh1{Ubar}; GInh_hpf1(r,prePnts(b+1)+BarPnts(b+1)*(i(b+1,a)-1):prePnts(b+1)+BarPnts(b+1)*i(b+1,a))] ;
        gAlt_Exc1{Ubar} = [gAlt_Exc1{Ubar}; GAlt_ExcInt_hpf1(r,prePnts(b+2)+BarPnts(b+2)*(i(b+2,a)-1):prePnts(b+2)+BarPnts(b+2)*i(b+2,a))] ;
        gAlt_Inh1{Ubar} = [gAlt_Inh1{Ubar}; GAlt_InhInt_hpf1(r,prePnts(b+2)+BarPnts(b+2)*(i(b+2,a)-1):prePnts(b+2)+BarPnts(b+2)*i(b+2,a))] ;
    
        gExc2{Ubar} = [gExc2{Ubar}; GExc_hpf2(r,prePnts(b)+BarPnts(b)*(i(b,a)-1):prePnts(b)+BarPnts(b)*i(b,a))] ;
        gInh2{Ubar} = [gInh2{Ubar}; GInh_hpf2(r,prePnts(b+1)+BarPnts(b+1)*(i(b+1,a)-1):prePnts(b+1)+BarPnts(b+1)*i(b+1,a))] ; 
        gAlt_Exc2{Ubar} = [gAlt_Exc2{Ubar}; GAlt_ExcInt_hpf2(r,prePnts(b+2)+BarPnts(b+2)*(i(b+2,a)-1):prePnts(b+2)+BarPnts(b+2)*i(b+2,a))] ;
        gAlt_Inh2{Ubar} = [gAlt_Inh2{Ubar}; GAlt_InhInt_hpf2(r,prePnts(b+2)+BarPnts(b+2)*(i(b+2,a)-1):prePnts(b+2)+BarPnts(b+2)*i(b+2,a))] ;    
    end
end

numBarRepeats = size(gExc1{1},1) ;

% % cut up and arrange data in array by bar direction and on/off responses
% r=0;
% [Ba,i] = sort(BarAngles,2) ; % sort each row of bar angle in accending order
% for b= 1:3:size(BarAngles,1) ; % on each trial
%     r=r+1 ;
%     for  a = 1:NumBars(1) ; % for each bar
%         gExc1{a}(r,:) = GExc_hpf1(r,prePnts(b)+BarPnts(b)*(i(b,a)-1):prePnts(b)+BarPnts(b)*i(b,a)) ;
%         gInh1{a}(r,:) = GInh_hpf1(r,prePnts(b+1)+BarPnts(b+1)*(i(b+1,a)-1):prePnts(b+1)+BarPnts(b+1)*i(b+1,a)) ;
%         gAlt_Exc1{a}(r,:) = GAlt_ExcInt_hpf1(r,prePnts(b+2)+BarPnts(b+2)*(i(b+2,a)-1):prePnts(b+2)+BarPnts(b+2)*i(b+2,a)) ;
%         gAlt_Inh1{a}(r,:) = GAlt_InhInt_hpf1(r,prePnts(b+2)+BarPnts(b+2)*(i(b+2,a)-1):prePnts(b+2)+BarPnts(b+2)*i(b+2,a)) ;
%     
%         gExc2{a}(r,:) = GExc_hpf2(r,prePnts(b)+BarPnts(b)*(i(b,a)-1):prePnts(b)+BarPnts(b)*i(b,a)) ;
%         gInh2{a}(r,:) = GInh_hpf2(r,prePnts(b+1)+BarPnts(b+1)*(i(b+1,a)-1):prePnts(b+1)+BarPnts(b+1)*i(b+1,a)) ;
%         gAlt_Exc2{a}(r,:) = GAlt_ExcInt_hpf2(r,prePnts(b+2)+BarPnts(b+2)*(i(b+2,a)-1):prePnts(b+2)+BarPnts(b+2)*i(b+2,a)) ;
%         gAlt_Inh2{a}(r,:) = GAlt_InhInt_hpf2(r,prePnts(b+2)+BarPnts(b+2)*(i(b+2,a)-1):prePnts(b+2)+BarPnts(b+2)*i(b+2,a)) ;    
%     end
% end

% get mean data
for a = 1:NumUniqueBars ;
    gExc_Mean1{a} = mean(gExc1{a}) ;
    gInh_Mean1{a} = mean(gInh1{a}) ;
    gAlt_Exc_Mean1{a} = mean(gAlt_Exc1{a}) ;
    gAlt_Inh_Mean1{a} = mean(gAlt_Inh1{a}) ;
    
    gExc_Mean2{a} = mean(gExc2{a}) ;
    gInh_Mean2{a} = mean(gInh2{a}) ;
    gAlt_Exc_Mean2{a} = mean(gAlt_Exc2{a}) ;
    gAlt_Inh_Mean2{a} = mean(gAlt_Inh2{a}) ;
end

% get residuals
residualOption = 0 ;
for a = 1:NumUniqueBars ;
    if residualOption == 0 ;
        gExc_Res1{a} =  gExc1{a} - repmat(gExc_Mean1{a},numBarRepeats,1) ;
        gInh_Res1{a} =  gInh1{a} - repmat(gInh_Mean1{a},numBarRepeats,1) ;
        gAlt_Exc_Res1{a} = gAlt_Exc1{a} - repmat(gAlt_Exc_Mean1{a},numBarRepeats,1) ;
        gAlt_Inh_Res1{a} = gAlt_Inh1{a} - repmat(gAlt_Inh_Mean1{a},numBarRepeats,1) ;

        gExc_Res2{a} =  gExc2{a} - repmat(gExc_Mean2{a},numBarRepeats,1) ;
        gInh_Res2{a} =  gInh2{a} - repmat(gInh_Mean2{a},numBarRepeats,1) ;
        gAlt_Exc_Res2{a} = gAlt_Exc2{a} - repmat(gAlt_Exc_Mean2{a},numBarRepeats,1) ;
        gAlt_Inh_Res2{a} = gAlt_Inh2{a} - repmat(gAlt_Inh_Mean2{a},numBarRepeats,1) ;
    
    elseif residualOption == 1 ;
        for b = 2:numBarRepeats-1 ;
            gExc_Res1{a}(b-1,:) =  gExc1{a}(b,:) - (gExc1{a}(b-1,:)+gExc1{a}(b+1,:))/2 ;
            gInh_Res1{a}(b-1,:) =  gInh1{a}(b,:) - (gInh1{a}(b-1,:)+gInh1{a}(b+1,:))/2 ;
            gAlt_Exc_Res1{a}(b-1,:) =  gAlt_Exc1{a}(b,:) - (gAlt_Exc1{a}(b-1,:)+gAlt_Exc1{a}(b+1,:))/2 ;
            gAlt_Inh_Res1{a}(b-1,:) =  gAlt_Inh1{a}(b,:) - (gAlt_Inh1{a}(b-1,:)+gAlt_Inh1{a}(b+1,:))/2 ;

            gExc_Res2{a}(b-1,:) =  gExc2{a}(b,:) - (gExc2{a}(b-1,:)+gExc2{a}(b+1,:))/2 ;
            gInh_Res2{a}(b-1,:) =  gInh2{a}(b,:) - (gInh2{a}(b-1,:)+gInh2{a}(b+1,:))/2 ;
            gAlt_Exc_Res2{a}(b-1,:) =  gAlt_Exc2{a}(b,:) - (gAlt_Exc2{a}(b-1,:)+gAlt_Exc2{a}(b+1,:))/2 ;
            gAlt_Inh_Res2{a}(b-1,:) =  gAlt_Inh2{a}(b,:) - (gAlt_Inh2{a}(b-1,:)+gAlt_Inh2{a}(b+1,:))/2 ;
        end
    end    
end
                
% g peaks for tunning curves 
for a= 1:NumUniqueBars ;
    gExc_max_Mean1(a) = mean(max(gExc1{a},[],2)) ;
    gInh_max_Mean1(a) = mean(max(gInh1{a},[],2)) ;
    gAlt_Exc_max_Mean1(a) = mean(max(gAlt_Exc1{a},[],2)) ;
    gAlt_Inh_max_Mean1(a) = mean(max(gAlt_Inh1{a},[],2)) ;
    
    gExc_max_Mean2(a) = mean(max(gExc2{a},[],2)) ;
    gInh_max_Mean2(a) = mean(max(gInh2{a},[],2)) ;
    gAlt_Exc_max_Mean2(a) = mean(max(gAlt_Exc2{a},[],2)) ;
    gAlt_Inh_max_Mean2(a) = mean(max(gAlt_Inh2{a},[],2)) ;
    
    gExc_max_std1(a) = std(max(gExc1{a},[],2)) ;
    gInh_max_std1(a) = std(max(gInh1{a},[],2)) ;
    gAlt_Exc_max_std1(a) = std(max(gAlt_Exc1{a},[],2)) ;
    gAlt_Inh_max_std1(a) = std(max(gAlt_Inh1{a},[],2)) ;
    
    gExc_max_std2(a) = std(max(gExc2{a},[],2)) ;
    gInh_max_std2(a) = std(max(gInh2{a},[],2)) ;
    gAlt_Exc_max_std2(a) = std(max(gAlt_Exc2{a},[],2)) ;
    gAlt_Inh_max_std2(a) = std(max(gAlt_Inh2{a},[],2)) ;
end
    
%time vectors for cc
time_cc = [SI(1)*([1:2*length(gAlt_Exc1{a})-1] - length(gAlt_Exc1{a}))] ;
   
% get cross correlations 
for a= 1:NumUniqueBars ;
    cc1(a,:) = xcov(gExc_Mean1{a},gInh_Mean1{a},'coef') ;
    ccAlt1(a,:) = xcov(gAlt_Exc_Mean1{a},gAlt_Inh_Mean1{a},'coef') ;
    
    cc2(a,:) = xcov(gExc_Mean2{a},gInh_Mean2{a},'coef') ;
    ccAlt2(a,:) = xcov(gAlt_Exc_Mean2{a},gAlt_Inh_Mean2{a},'coef') ;
       
    for b = 1:numBarRepeats-residualOption*2;    
        ccRes1{a}(b,:) = xcov(gExc_Res1{a}(b,:),gInh_Res1{a}(b,:),'coef') ; % cc of res inh1/exc1
        ccAltRes1{a}(b,:) = xcov(gAlt_Exc_Res1{a}(b,:),gAlt_Inh_Res1{a}(b,:),'coef') ;

        ccRes2{a}(b,:) = xcov(gExc_Res2{a}(b,:),gInh_Res2{a}(b,:),'coef') ; % cc of res
        ccAltRes2{a}(b,:) = xcov(gAlt_Exc_Res2{a}(b,:),gAlt_Inh_Res2{a}(b,:),'coef') ;
 
        ccExcRes{a}(b,:) = xcov(gExc_Res1{a}(b,:),gExc_Res2{a}(b,:),'coef') ; % cc of res exc1/exc2
        ccAltExcRes{a}(b,:) = xcov(gAlt_Exc_Res1{a}(b,:),gAlt_Exc_Res2{a}(b,:),'coef') ; 
        
        ccInhRes{a}(b,:) = xcov(gInh_Res1{a}(b,:),gInh_Res2{a}(b,:),'coef') ; % cc of res inh1/inh2 
        ccAltInhRes{a}(b,:) = xcov(gAlt_Inh_Res1{a}(b,:),gAlt_Inh_Res2{a}(b,:),'coef') ; 
        
        % estimate of converging correlations if simultaneous recording was not possible
        ccE1E2{a}(b,:) = xcov(gAlt_Exc_Res1{a}(b,:),gAlt_Exc_Res2{a}(b,:)) ; % cc of res exc/exc
        ccI1I2{a}(b,:) = xcov(gAlt_Inh_Res1{a}(b,:),gAlt_Inh_Res2{a}(b,:)) ; % cc of res inh/inh
        ccE1I2{a}(b,:) = xcov(gAlt_Exc_Res1{a}(b,:),gAlt_Inh_Res2{a}(b,:)) ; % cc of res exc/inh
        ccE2I1{a}(b,:) = xcov(gAlt_Exc_Res2{a}(b,:),gAlt_Inh_Res1{a}(b,:)) ; % cc of res exc/inh
        
        acE1{a}(b,:) = xcov(gAlt_Exc_Res1{a}(b,:)) ; % autocorr
        acE2{a}(b,:) = xcov(gAlt_Exc_Res2{a}(b,:)) ;
        acI1{a}(b,:) = xcov(gAlt_Inh_Res1{a}(b,:)) ;
        acI2{a}(b,:) = xcov(gAlt_Inh_Res2{a}(b,:)) ;        
        
    end
    
    % mean cross corr
    ccResMean1(a,:) = mean(ccRes1{a}) ;
    ccAltResMean1(a,:) = mean(ccAltRes1{a}) ;

    ccResMean2(a,:) = mean(ccRes2{a}) ;
    ccAltResMean2(a,:) = mean(ccAltRes2{a}) ;    
    
    ccExcResMean(a,:) = mean(ccExcRes{a}) ;
    ccAltExcResMean(a,:) = mean(ccAltExcRes{a}) ;
    
    ccInhResMean(a,:) = mean(ccInhRes{a}) ;
    ccAltInhResMean(a,:) = mean(ccAltInhRes{a}) ;
    
    ccE1E2Mean(a,:) = mean(ccE1E2{a}) ;
    ccI1I2Mean(a,:) = mean(ccI1I2{a}) ;
    ccEIMean(a,:) = mean([ccE1I2{a};ccE2I1{a}]) ;
    
    acE1Mean(a,:) = mean(acE1{a}) ;
    acE2Mean(a,:) = mean(acE2{a}) ;
    acI1Mean(a,:) = mean(acI1{a}) ;
    acI2Mean(a,:) = mean(acI2{a}) ;

    % standard error of mean (sem) of cross correlations 
    ccRes_SEM1(a,:) = std(ccRes1{a})/sqrt(numBarRepeats) ;
    ccAltRes_SEM1(a,:) = std(ccAltRes1{a})/sqrt(numBarRepeats) ;
   
    ccRes_SEM2(a,:) = std(ccRes2{a})/sqrt(numBarRepeats) ;
    ccAltRes_SEM2(a,:) = std(ccAltRes2{a})/sqrt(numBarRepeats) ;
    
    % cc peaks
    ccAltResPeak1(a) = CCpeakFinder(ccAltResMean1(a,:)) ;
    ccAltResPeak2(a) = CCpeakFinder(ccAltResMean2(a,:)) ;
    
    ccE1E2Mean_peak(a) = CCpeakFinder(ccE1E2Mean(a,:)) ;
    ccI1I2Mean_peak(a) = CCpeakFinder(ccI1I2Mean(a,:)) ; 
    ccEIMean_peak(a) = CCpeakFinder(ccEIMean(a,:)) ;
    
    acE1Mean_peak(a) = CCpeakFinder(acE1Mean(a,:)) ;
    acE2Mean_peak(a) = CCpeakFinder(acE2Mean(a,:)) ;
    acI1Mean_peak(a) = CCpeakFinder(acI1Mean(a,:)) ;
    acI2Mean_peak(a) = CCpeakFinder(acI2Mean(a,:)) ;
    
end

% estimate converging correlations from intercorrelations
Estimate = ccEIMean_peak./sqrt(abs(ccE1E2Mean_peak.*ccI1I2Mean_peak)) ;
lowerBound1 = ccEIMean_peak./sqrt(acE1Mean_peak.*acI1Mean_peak) ;
lowerBound2 = ccEIMean_peak./sqrt(acE2Mean_peak.*acI2Mean_peak) ;

upperBound1 = (ccEIMean_peak + min([(acE1Mean_peak-ccE1E2Mean_peak);(acI1Mean_peak-ccI1I2Mean_peak)],[],1))./sqrt(acE1Mean_peak.*acI1Mean_peak)  ;
upperBound2 = (ccEIMean_peak + min([(acE2Mean_peak-ccE1E2Mean_peak);(acI2Mean_peak-ccI1I2Mean_peak)],[],1))./sqrt(acE2Mean_peak.*acI2Mean_peak) ;


% % % figures
% figure
% plot(time,Monitor)
% % 
% % 
% figure % conductance check
% for a=1:numTrials ; % on each trial
%     subplot(4,1,1)
%     plot(time,GExc_hpf1(a,:))
%     hold on
%     plot(time,GInh_hpf1(a,:),'r')
%     hold off
%     
%     subplot(4,1,2)
%     plot(time,GExc_hpf2(a,:))
%     hold on
%     plot(time,GInh_hpf2(a,:),'r')
%     hold off
%     
%     subplot(4,1,3)
%     plot(time,GAlt_ExcInt_hpf1(a,:))
%     hold on
%     plot(time,GAlt_InhInt_hpf1(a,:),'r')
%     hold off
%     
%     subplot(4,1,4)
%     plot(time,GAlt_ExcInt_hpf2(a,:))
%     hold on
%     plot(time,GAlt_InhInt_hpf2(a,:),'r')
%     hold off
%     
%     text(0,.9,num2str(epochs((a-1)*3+1)),'units','norm')
%     pause
% end
% 
% 
figure
for a=1:NumUniqueBars ;
    subplot(ceil(NumUniqueBars(1)/2),2,a)
    plot(time_cc,ccAltResMean1(a,:),'k')
    hold on
    plot(time_cc,ccResMean1(a,:),'g')
end

figure
for a=1:NumUniqueBars ;
    subplot(ceil(NumUniqueBars(1)/2),2,a)
    plot(time_cc,ccAltResMean2(a,:),'k')
    hold on
    plot(time_cc,ccResMean2(a,:),'g')
end

figure
for a=1:NumUniqueBars ;
    subplot(ceil(NumUniqueBars(1)/2),2,a)
    plot(time_cc,ccAltExcResMean(a,:),'b')
    hold on
    plot(time_cc,ccAltInhResMean(a,:),'r')
end

figure
for a=1:NumUniqueBars ;
    subplot(ceil(NumUniqueBars(1)/2),2,a)
    plot(time_cc,ccExcResMean(a,:),'b')
    hold on
    plot(time_cc,ccInhResMean(a,:),'r')
end
% 
% figure
% plot(lowerBound1)
% hold on
% plot(upperBound1)
% plot(Estimate,'b--')
% plot(ccAltResPeak1,'k')
% plot(ccAltResPeak2,'k--')

% For igor

% prep G for dynamic clamp

maxExc1 = max(max(GAlt_ExcInt_hpf1(:,prePnts(1):end-postPnts(1)))) ;
maxExc2 = max(max(GAlt_ExcInt_hpf2(:,prePnts(1):end-postPnts(1)))) ;

maxInh1 = max(max(GAlt_InhInt_hpf1(:,prePnts(1):end-postPnts(1)))) ;
maxInh2 = max(max(GAlt_InhInt_hpf2(:,prePnts(1):end-postPnts(1)))) ;

for b=1:numBarRepeats ;
    ExcG1 = [] ;
    InhG1 = [] ;

    ExcG2 = [] ;
    InhG2 = [] ;

    for a=1:NumUniqueBars ;
        ExcG1 = [ExcG1,gAlt_Exc1{a}(b,:)] ;
        InhG1 = [InhG1,gAlt_Inh1{a}(b,:)] ;
        
        ExcG2 = [ExcG2,gAlt_Exc2{a}(b,:)] ;
        InhG2 = [InhG2,gAlt_Inh2{a}(b,:)] ;
    end
    
   % normalized for dynamic clamp
    identifier = ['ExcG1t',num2str(b),'c',num2str(A)] ;
    ForIgor.(identifier) = ExcG1/maxExc1 ;
    
    identifier = ['InhG1t',num2str(b),'c',num2str(A)] ;
    ForIgor.(identifier) = InhG1/maxInh1 ;
       
    identifier = ['ExcG2t',num2str(b),'c',num2str(A)] ;
    ForIgor.(identifier) = ExcG2/maxExc2 ;
       
    identifier = ['InhG2t',num2str(b),'c',num2str(A)] ;
    ForIgor.(identifier) = InhG2/maxInh2 ;


%     % conductances not normalized
%     identifier = ['ExcG1t',num2str(b),'c',num2str(A)] ;
%     ForIgor.(identifier) = ExcG1 ;
%     
%     identifier = ['InhG1t',num2str(b),'c',num2str(A)] ;
%     ForIgor.(identifier) = InhG1 ;
%        
%     identifier = ['ExcG2t',num2str(b),'c',num2str(A)] ;
%     ForIgor.(identifier) = ExcG2 ;
%        
%     identifier = ['InhG2t',num2str(b),'c',num2str(A)] ;
%     ForIgor.(identifier) = InhG2 ;
end


% % cross correlations
% identifier = ['timecc',num2str(A)] ;
% ForIgor.(identifier) = time_cc ;
% 
% for a=1:NumUniqueBars ;
%     identifier = ['ccAltResMean1b',num2str(a),'c',num2str(A)] ;
%     ForIgor.(identifier) = ccAltResMean1(a,:) ;
% 
%     identifier = ['ccResMean1b',num2str(a),'c',num2str(A)] ;
%     ForIgor.(identifier) = ccResMean1(a,:) ;
% 
%     identifier = ['ccAltResMean2b',num2str(a),'c',num2str(A)] ;
%     ForIgor.(identifier) = ccAltResMean2(a,:) ;
% 
%     identifier = ['ccResMean2b',num2str(a),'c',num2str(A)] ;
%     ForIgor.(identifier) = ccResMean2(a,:) ;
% 
%     identifier = ['ccAltExcResMeanb',num2str(a),'c',num2str(A)] ;
%     ForIgor.(identifier) = ccAltExcResMean(a,:) ;
% 
%     identifier = ['ccAltInhResMeanb',num2str(a),'c',num2str(A)] ;
%     ForIgor.(identifier) = ccAltInhResMean(a,:) ;
% 
%     identifier = ['ccExcResMeanb',num2str(a),'c',num2str(A)] ;
%     ForIgor.(identifier) = ccExcResMean(a,:) ;
%     
%     identifier = ['ccInhResMeanb',num2str(a),'c',num2str(A)] ;
%     ForIgor.(identifier) = ccInhResMean(a,:) ;
% end






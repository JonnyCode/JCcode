function ForIgor = ValtAnalyzerDS2(Input,Parameters,id,id2,A) ; 
% this function will analyze alternate voltage exp from ds cells presented
% a moving bar

% JC 8/7/09 adapted from earlier version
% editted 7/20/11 to withhold trial from mean used to calculate the residual for that trial

residualOption = Parameters.residualOption ;

% get data
[fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/mouse/',Input(A).cellname]);

epochs = str2num(Input(A).(id)) ;
round = 0 ;
for a = 1:3:length(epochs) ; % for each spike epoch
    round = round +1 ;

    [dataExc(round,:), error] = ITCReadEpoch(epochs(a), 0, fp) ;    %#ok<*AGROW> % get data

    [dataInh(round,:), error] = ITCReadEpoch(epochs(a+1), 0, fp) ;    % get data
    
    [dataAltV(round,:), error] = ITCReadEpoch(epochs(a+2), 0, fp) ;    % get data
    
    [voltageCommand(round,:), error] = ITCReadEpochStm(epochs(a+2), 0,fp); % get voltage command
    
    [SI(round), error] = ITCGetSamplingInterval(epochs(a+2), fp); % get sampling interval
    SI(round) = SI(round) * 1e-6; % Sampling interval in sec
end

voltageCommand = voltageCommand(:,1:length(dataExc)) ;

if Input(A).ITC18flag == 1 ;
    SI = SI*1.25 ;
end

SpatialStimParams = hdf5load(['~/Data/mouse/',Input(A).cellname,'_spatial.h5']) ; % load spatial stim params

for a = 1:length(epochs) ;
    [Monitor(a,:), error] = ITCReadEpoch(epochs(a), 1, fp) ;    % get monitor data
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

    stimPnts(a) = floor(Struct.spatial_stimpts/(frameRate*SI(1))) ;
    prePnts(a) = floor(Struct.spatial_prepts/(frameRate*SI(1))) ;
    postPnts(a) = floor(Struct.spatial_postpts/(frameRate*SI(1))) ;
 
    NumBars(a) = length(BarAngles(a,:)) ; % number of bars shown per trial
    %BarPnts(a) = floor((length(dataExc)-prePnts(a)-postPnts(a))/NumBars(a)) ; % number of points the bar is presented + interbar points PREVIOUS CODE USE BEFORE 7/15/11 PROBLEM WITH LENGTH OF DATA IGOR IS RECORDING
    BarPnts(a) = floor(stimPnts(a)/NumBars(a)) ;
    OnPnts(a) = floor((BarWidth(a)/BarSpeed(a))/(frameRate*SI(1))) ; % number of points before end of bar appears triggering off response
end

cyclepnts = 100 ; % number of sample points in a cycle, pnts between leaving hold1 and returning (Also gets rid of first cycle)    
FirstAltPnt = (cyclepnts/2)+1 ; % first sample point you want to plot after begining of step from alternation 
%FirstAltPnt = 30 ;
LastAltPnt = FirstAltPnt ;  % last "                                                                    "

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

% low pass filter and remove electrical crap
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

offsetExc = min(min([GExc_hpf(:,prePnts:prePnts+stimPnts),GAlt_ExcInt_hpf(:,prePnts:prePnts+stimPnts)])) ;
offsetInh = min(min([GInh_hpf(:,prePnts:prePnts+stimPnts),GAlt_InhInt_hpf(:,prePnts:prePnts+stimPnts)])) ;
% ofset G (assumes all g have same mean and 1 min) 
for a = 1:size(dataExc,1) ; % for each trial
    GExc_hpf(a,:) = GExc_hpf(a,:) - offsetExc ; % offsets
    GInh_hpf(a,:) = GInh_hpf(a,:) - offsetInh ; 
    GAlt_ExcInt_hpf(a,:) = GAlt_ExcInt_hpf(a,:) - offsetExc ; 
    GAlt_InhInt_hpf(a,:) = GAlt_InhInt_hpf(a,:) - offsetInh ;
end

% to save memmory clear unused variables
clearvars -except GExc_hpf GInh_hpf GAlt_ExcInt_hpf GAlt_InhInt_hpf...
    prePnts BarPnts OnPnts BarAngles NumBars SI Monitor A Input fp...
    frameRate SpatialStimParams id id2 time residualOption

% cut up and arrange data in array by bar direction and on/off responses
r=0;
[Ba,i] = sort(BarAngles,2) ; % sort each row of bar angle in accending order
for b= 1:3:size(BarAngles,1) ; % on each trial
    r=r+1 ;
    for  a = 1:NumBars(1) ; % for each bar
        
        gExc{a}(r,:) = GExc_hpf(r,prePnts(b)+BarPnts(b)*(i(b,a)-1):prePnts(b)+BarPnts(b)*i(b,a)) ;
        gInh{a}(r,:) = GInh_hpf(r,prePnts(b+1)+BarPnts(b+1)*(i(b+1,a)-1):prePnts(b+1)+BarPnts(b+1)*i(b+1,a)) ;
        gAlt_Exc{a}(r,:) = GAlt_ExcInt_hpf(r,prePnts(b+2)+BarPnts(b+2)*(i(b+2,a)-1):prePnts(b+2)+BarPnts(b+2)*i(b+2,a)) ;
        gAlt_Inh{a}(r,:) = GAlt_InhInt_hpf(r,prePnts(b+2)+BarPnts(b+2)*(i(b+2,a)-1):prePnts(b+2)+BarPnts(b+2)*i(b+2,a)) ;
    
        gExc_ON{a}(r,:) = GExc_hpf(r,prePnts(b)+BarPnts(b)*(i(b,a)-1):prePnts(b)+BarPnts(b)*(i(b,a)-1)+OnPnts(b)) ;
        gInh_ON{a}(r,:) = GInh_hpf(r,prePnts(b+1)+BarPnts(b+1)*(i(b+1,a)-1):prePnts(b+1)+BarPnts(b+1)*(i(b+1,a)-1)+OnPnts(b+1)) ;
        gAlt_Exc_ON{a}(r,:) = GAlt_ExcInt_hpf(r,prePnts(b+2)+BarPnts(b+2)*(i(b+2,a)-1):prePnts(b+2)+BarPnts(b+2)*(i(b+2,a)-1)+OnPnts(b+2)) ;
        gAlt_Inh_ON{a}(r,:) = GAlt_InhInt_hpf(r,prePnts(b+2)+BarPnts(b+2)*(i(b+2,a)-1):prePnts(b+2)+BarPnts(b+2)*(i(b+2,a)-1)+OnPnts(b+2)) ;
    
        gExc_OFF{a}(r,:) = GExc_hpf(r,prePnts(b)+BarPnts(b)*(i(b,a)-1)+OnPnts(b):prePnts(b)+BarPnts(b)*(i(b,a))) ;
        gInh_OFF{a}(r,:) = GInh_hpf(r,prePnts(b+1)+BarPnts(b+1)*(i(b+1,a)-1)+OnPnts(b+1):prePnts(b+1)+BarPnts(b+1)*(i(b+1,a))) ;
        gAlt_Exc_OFF{a}(r,:) = GAlt_ExcInt_hpf(r,prePnts(b+2)+BarPnts(b+2)*(i(b+2,a)-1)+OnPnts(b+2):prePnts(b+2)+BarPnts(b+2)*(i(b+2,a))) ;
        gAlt_Inh_OFF{a}(r,:) = GAlt_InhInt_hpf(r,prePnts(b+2)+BarPnts(b+2)*(i(b+2,a)-1)+OnPnts(b+2):prePnts(b+2)+BarPnts(b+2)*(i(b+2,a))) ;
    
    end
end

% get mean, meanamp, and mean time peak  data
for a = 1:NumBars(1) ;
    gExc_Mean{a} = mean(gExc{a}) ;
    gInh_Mean{a} = mean(gInh{a}) ;
    gAlt_Exc_Mean{a} = mean(gAlt_Exc{a}) ;
    gAlt_Inh_Mean{a} = mean(gAlt_Inh{a}) ;
    
    gExc_ON_Mean{a} = mean(gExc_ON{a}) ;
    gInh_ON_Mean{a} = mean(gInh_ON{a}) ;
    gAlt_Exc_ON_Mean{a} = mean(gAlt_Exc_ON{a}) ;
    gAlt_Inh_ON_Mean{a} = mean(gAlt_Inh_ON{a}) ;
    
    gExc_OFF_Mean{a} = mean(gExc_OFF{a}) ;
    gInh_OFF_Mean{a} = mean(gInh_OFF{a}) ;
    gAlt_Exc_OFF_Mean{a} = mean(gAlt_Exc_OFF{a}) ;
    gAlt_Inh_OFF_Mean{a} = mean(gAlt_Inh_OFF{a}) ;

    gExc_ON_meanamp{a} = mean(gExc_ON{a},2) ; % mean amp
    gInh_ON_meanamp{a} = mean(gInh_ON{a},2) ;
    
    for b=1:size(gExc_ON{a},1);
        gExc_ON_timemaxamp{a}(b) = time(gExc_ON{a}(b,:) == max(gExc_ON{a}(b,:))) ; % time of max amp
        gInh_ON_timemaxamp{a}(b) = time(gInh_ON{a}(b,:) == max(gInh_ON{a}(b,:))) ;
    end
    
end

% get residuals
if residualOption == 0 ; % added residual option 8/4/10
    % standard
    for a= 1:NumBars(1) ;  
       for b=1:size(gExc{a},1) ; % for each possible residual trial
            gExc_Res{a}(b,:) =  gExc{a}(b,:) - mean(gExc{a})+(gExc{a}(b,:)/size(gExc{a},1)) ;
            gInh_Res{a}(b,:) =  gInh{a}(b,:) - mean(gInh{a})+(gInh{a}(b,:)/size(gExc{a},1)) ;
            gAlt_Exc_Res{a}(b,:) = gAlt_Exc{a}(b,:) - mean(gAlt_Exc{a})+(gAlt_Exc{a}(b,:)/size(gExc{a},1)) ;
            gAlt_Inh_Res{a}(b,:) = gAlt_Inh{a}(b,:) - mean(gAlt_Inh{a})+(gAlt_Inh{a}(b,:)/size(gExc{a},1)) ;

            gExc_ON_Res{a}(b,:) =  gExc_ON{a}(b,:) - mean(gExc_ON{a})+(gExc_ON{a}(b,:)/size(gExc{a},1)) ;
            gInh_ON_Res{a}(b,:) =  gInh_ON{a}(b,:) - mean(gInh_ON{a})+(gInh_ON{a}(b,:)/size(gExc{a},1)) ;
            gAlt_Exc_ON_Res{a}(b,:) = gAlt_Exc_ON{a}(b,:) - mean(gAlt_Exc_ON{a})+(gAlt_Exc_ON{a}(b,:)/size(gExc{a},1)) ;
            gAlt_Inh_ON_Res{a}(b,:) = gAlt_Inh_ON{a}(b,:) - mean(gAlt_Inh_ON{a})+(gAlt_Inh_ON{a}(b,:)/size(gExc{a},1)) ;

            gExc_OFF_Res{a}(b,:) =  gExc_OFF{a}(b,:) - mean(gExc_OFF{a})+(gExc_OFF{a}(b,:)/size(gExc{a},1)) ;
            gInh_OFF_Res{a}(b,:) =  gInh_OFF{a}(b,:) - mean(gInh_OFF{a})+(gInh_OFF{a}(b,:)/size(gExc{a},1)) ;
            gAlt_Exc_OFF_Res{a}(b,:) = gAlt_Exc_OFF{a}(b,:) - mean(gAlt_Exc_OFF{a})+(gAlt_Exc_OFF{a}(b,:)/size(gExc{a},1)) ;
            gAlt_Inh_OFF_Res{a}(b,:) = gAlt_Inh_OFF{a}(b,:) - mean(gAlt_Inh_OFF{a})+(gAlt_Inh_OFF{a}(b,:)/size(gExc{a},1)) ;
         
        end
    end

elseif residualOption == 1 ;
    % nearest neighbor residuals
    for a= 1:NumBars(1) ;
        for b=2:size(gExc{a},1)-1 ; % for each possible residual trial
            gExc_Res{a}(b-1,:) =  gExc{a}(b,:) - (gExc{a}(b-1,:)+gExc{a}(b+1,:))/2 ;
            gInh_Res{a}(b-1,:) =  gInh{a}(b,:) - (gInh{a}(b-1,:)+gInh{a}(b+1,:))/2 ;
            gAlt_Exc_Res{a}(b-1,:) = gAlt_Exc{a}(b,:) - (gAlt_Exc{a}(b-1,:)+gAlt_Exc{a}(b+1,:))/2 ;
            gAlt_Inh_Res{a}(b-1,:) = gAlt_Inh{a}(b,:) - (gAlt_Inh{a}(b-1,:)+gAlt_Inh{a}(b+1,:))/2 ;

            gExc_ON_Res{a}(b-1,:) =  gExc_ON{a}(b,:) - (gExc_ON{a}(b-1,:)+gExc_ON{a}(b+1,:))/2 ;
            gInh_ON_Res{a}(b-1,:) =  gInh_ON{a}(b,:) - (gInh_ON{a}(b-1,:)+gInh_ON{a}(b+1,:))/2 ;
            gAlt_Exc_ON_Res{a}(b-1,:) = gAlt_Exc_ON{a}(b,:) - (gAlt_Exc_ON{a}(b-1,:)+gAlt_Exc_ON{a}(b+1,:))/2 ;
            gAlt_Inh_ON_Res{a}(b-1,:) = gAlt_Inh_ON{a}(b,:) - (gAlt_Inh_ON{a}(b-1,:)+gAlt_Inh_ON{a}(b+1,:))/2 ;

            gExc_OFF_Res{a}(b-1,:) =  gExc_OFF{a}(b,:) - (gExc_OFF{a}(b-1,:)+gExc_OFF{a}(b+1,:))/2 ;
            gInh_OFF_Res{a}(b-1,:) =  gInh_OFF{a}(b,:) - (gInh_OFF{a}(b-1,:)+gInh_OFF{a}(b+1,:))/2 ;
            gAlt_Exc_OFF_Res{a}(b-1,:) = gAlt_Exc_OFF{a}(b,:) - (gAlt_Exc_OFF{a}(b-1,:)+gAlt_Exc_OFF{a}(b+1,:))/2 ;
            gAlt_Inh_OFF_Res{a}(b-1,:) = gAlt_Inh_OFF{a}(b,:) - (gAlt_Inh_OFF{a}(b-1,:)+gAlt_Inh_OFF{a}(b+1,:))/2 ;
         
%             gExc_ON_meanamp_Res{a}(b-1,:) =  gExc_ON_meanamp{a}(b,:) - (gExc_ON_meanamp{a}(b-1,:)+gExc_ON_meanamp{a}(b+1,:))/2 ;
%             gInh_ON_meanamp_Res{a}(b-1,:) =  gInh_ON_meanamp{a}(b,:) - (gInh_ON_meanamp{a}(b-1,:)+gInh_ON_meanamp{a}(b+1,:))/2 ;
        end
        
    end
end
    
    
% g peaks for tunning curves 
for a= 1:NumBars(1) ;
    gExc_max_Mean(a) = mean(max(gExc{a},[],2)) ;
    gInh_max_Mean(a) = mean(max(gInh{a},[],2)) ;
    gAlt_Exc_max_Mean(a) = mean(max(gAlt_Exc{a},[],2)) ;
    gAlt_Inh_max_Mean(a) = mean(max(gAlt_Inh{a},[],2)) ;
    
    gExc_ON_max_Mean(a) = mean(max(gExc_ON{a},[],2)) ;
    gInh_ON_max_Mean(a) = mean(max(gInh_ON{a},[],2)) ;
    gAlt_Exc_ON_max_Mean(a) = mean(max(gAlt_Exc_ON{a},[],2)) ;
    gAlt_Inh_ON_max_Mean(a) = mean(max(gAlt_Inh_ON{a},[],2)) ;
    
    gExc_OFF_max_Mean(a) = mean(max(gExc_OFF{a},[],2)) ;
    gInh_OFF_max_Mean(a) = mean(max(gInh_OFF{a},[],2)) ;
    gAlt_Exc_OFF_max_Mean(a) = mean(max(gAlt_Exc_OFF{a},[],2)) ;
    gAlt_Inh_OFF_max_Mean(a) = mean(max(gAlt_Inh_OFF{a},[],2)) ;
    
    gExc_max_std(a) = std(max(gExc{a},[],2)) ;
    gInh_max_std(a) = std(max(gInh{a},[],2)) ;
    gAlt_Exc_max_std(a) = std(max(gAlt_Exc{a},[],2)) ;
    gAlt_Inh_max_std(a) = std(max(gAlt_Inh{a},[],2)) ;
    
    gExc_ON_max_std(a) = std(max(gExc_ON{a},[],2)) ;
    gInh_ON_max_std(a) = std(max(gInh_ON{a},[],2)) ;
    gAlt_Exc_ON_max_std(a) = std(max(gAlt_Exc_ON{a},[],2)) ;
    gAlt_Inh_ON_max_std(a) = std(max(gAlt_Inh_ON{a},[],2)) ;
    
    gExc_OFF_max_std(a) = std(max(gExc_OFF{a},[],2)) ;
    gInh_OFF_max_std(a) = std(max(gInh_OFF{a},[],2)) ;
    gAlt_Exc_OFF_max_std(a) = std(max(gAlt_Exc_OFF{a},[],2)) ;
    gAlt_Inh_OFF_max_std(a) = std(max(gAlt_Inh_OFF{a},[],2)) ;
end
    
%time vectors for cc
time_cc = [SI(1)*([1:2*length(gAlt_Exc{a})-1] - length(gAlt_Exc{a}))] ;
time_ccON = [SI(1)*([1:2*length(gAlt_Exc_ON{a})-1] - length(gAlt_Exc_ON{a}))] ; 
time_ccOFF = [SI(1)*([1:2*length(gAlt_Exc_OFF{a})-1] - length(gAlt_Exc_OFF{a}))] ;    

% get cross correlations 
for a= 1:NumBars(1) ;
    cc(a,:) = xcov(gExc_Mean{a},gInh_Mean{a},'coef') ;
    ccAlt(a,:) = xcov(gAlt_Exc_Mean{a},gAlt_Inh_Mean{a},'coef') ;
    
    cc_ON(a,:) = xcov(gExc_ON_Mean{a},gInh_ON_Mean{a},'coef') ;
    ccAlt_ON(a,:) = xcov(gAlt_Exc_ON_Mean{a},gAlt_Inh_ON_Mean{a},'coef') ; 
    
    cc_OFF(a,:) = xcov(gExc_OFF_Mean{a},gInh_OFF_Mean{a},'coef') ;
    ccAlt_OFF(a,:) = xcov(gAlt_Exc_OFF_Mean{a},gAlt_Inh_OFF_Mean{a},'coef') ;
    
    gAlt_Inh_ResShuff{a} = circshift(gAlt_Inh_Res{a},[1,0]) ;
    gAlt_Inh_ON_ResShuff{a} = circshift(gAlt_Inh_ON_Res{a},[1,0]) ;
    gAlt_Inh_OFF_ResShuff{a} = circshift(gAlt_Inh_OFF_Res{a},[1,0]) ;
    
    for b = 1:size(gExc_Res{a},1) ;
    
        ccAltInd{a}(b,:) = xcov(gAlt_Exc{a}(b,:),gAlt_Inh{a}(b,:),'coef') ; % cc of the ind alt (these may be different than the cc of the mean g above)
        ccRes{a}(b,:) = xcov(gExc_Res{a}(b,:),gInh_Res{a}(b,:),'coef') ; % cc of res
        ccAltRes{a}(b,:) = xcov(gAlt_Exc_Res{a}(b,:),gAlt_Inh_Res{a}(b,:),'coef') ;
        ccAltResShuff{a}(b,:) = xcov(gAlt_Exc_Res{a}(b,:),gAlt_Inh_ResShuff{a}(b,:),'coef') ; % cc of res shuffled

        ccAltInd_ON{a}(b,:) = xcov(gAlt_Exc_ON{a}(b,:),gAlt_Inh_ON{a}(b,:),'coef') ;
        ccRes_ON{a}(b,:) = xcov(gExc_ON_Res{a}(b,:),gInh_ON_Res{a}(b,:),'coef') ;
        ccAltRes_ON{a}(b,:) = xcov(gAlt_Exc_ON_Res{a}(b,:),gAlt_Inh_ON_Res{a}(b,:),'coef') ;
        ccAltResShuff_ON{a}(b,:) = xcov(gAlt_Exc_ON_Res{a}(b,:),gAlt_Inh_ON_ResShuff{a}(b,:),'coef') ;

        ccAltInd_OFF{a}(b,:) = xcov(gAlt_Exc_OFF{a}(b,:),gAlt_Inh_OFF{a}(b,:),'coef') ;
        ccRes_OFF{a}(b,:) = xcov(gExc_OFF_Res{a}(b,:),gInh_OFF_Res{a}(b,:),'coef') ;
        ccAltRes_OFF{a}(b,:) = xcov(gAlt_Exc_OFF_Res{a}(b,:),gAlt_Inh_OFF_Res{a}(b,:),'coef') ;
        ccAltResShuff_OFF{a}(b,:) = xcov(gAlt_Exc_OFF_Res{a}(b,:),gAlt_Inh_OFF_ResShuff{a}(b,:),'coef') ;
    end
    % mean cross corr 
    ccAltIndMean(a,:) = mean(ccAltInd{a}) ;
    ccResMean(a,:) = mean(ccRes{a}) ;
    ccAltResMean(a,:) = mean(ccAltRes{a}) ;
    ccAltResShuffMean(a,:) = mean(ccAltResShuff{a}) ;
    
    ccAltInd_ON_Mean(a,:) = mean(ccAltInd_ON{a}) ;
    ccRes_ON_Mean(a,:) = mean(ccRes_ON{a}) ;
    ccAltRes_ON_Mean(a,:) = mean(ccAltRes_ON{a}) ;
    ccAltResShuff_ON_Mean(a,:) = mean(ccAltResShuff_ON{a}) ;
    
    ccAltInd_OFF_Mean(a,:) = mean(ccAltInd_OFF{a}) ;
    ccRes_OFF_Mean(a,:) = mean(ccRes_OFF{a}) ;
    ccAltRes_OFF_Mean(a,:) = mean(ccAltRes_OFF{a}) ;
    ccAltResShuff_OFF_Mean(a,:) = mean(ccAltResShuff_OFF{a}) ;
    
    % standard error of mean (sem) of cross correlations 
    ccRes_SEM(a,:) = std(ccRes{a})/sqrt(size(GExc_hpf,1)) ;
    ccAltRes_SEM(a,:) = std(ccAltRes{a})/sqrt(size(GExc_hpf,1)) ;
    ccAltResShuff_SEM(a,:) = std(ccAltResShuff{a})/sqrt(size(GExc_hpf,1)) ;    

    % cc peaks
    [ccRes_Peakabs,ccRes_Peaki] = max(abs(ccResMean(a,time_cc<.2 & time_cc>-.2))) ;
    ccRes_Peak(a) = ccResMean(a,ccRes_Peaki+find(time_cc>-.2,1)) ;

    [ccAltRes_Peakabs,ccAltRes_Peaki] = max(abs(ccAltResMean(a,time_cc<.2 & time_cc>-.2))) ;
    ccAltRes_Peak(a) = ccAltResMean(a,ccAltRes_Peaki+find(time_cc>-.2,1)) ;
    
%     % corr coef of meanamp
%     [cCoef,p] = corrcoef(gExc_ON_meanamp_Res{a},gInh_ON_meanamp_Res{a}) ;
%     corrCoef_meanamp(a) = cCoef(1,2) ;
%     p_meanamp(a) = p(1,2) ;
end
peakVector = [ccAltRes_Peak,nan(1,8-length(ccAltRes_Peak))] ; % make vector to align peaks
peakPlace = 4 ; % indicy of peak to which you should align
peakPlaceShift = peakPlace - find(ccAltRes_Peak == max(ccAltRes_Peak)) ;
peakVector = circshift(peakVector,[0,peakPlaceShift]) ;

% estimating covariance and variance of diveging noise assuming simple model
for a= 1:NumBars(1) ;
    for b = 1:size(gExc_Res{a},1) ;
        ccAltRes_forEstimate{a}(b,:) = xcov(gAlt_Exc_Res{a}(b,:),gAlt_Inh_Res{a}(b,:),'biased') ;
        ccAltRes_excAc_forEstimate{a}(b,:) = xcov(gAlt_Exc_Res{a}(b,:),'biased') ;
        ccAltRes_inhAc_forEstimate{a}(b,:) = xcov(gAlt_Inh_Res{a}(b,:),'biased') ;
    end
    
    ccAltRes_forEstimateMean(a,:) =  mean(ccAltRes_forEstimate{a}) ;
    [temp_Peakabs,temp_Peaki] = max(abs(ccAltRes_forEstimateMean(a,time_cc<.2 & time_cc>-.2))) ;
    ccAltRes_forEstimate_Peak(a) = ccAltRes_forEstimateMean(a,temp_Peaki+find(time_cc>-.2,1)) ; % covariance peak
    
    acAltExcRes_forEstimate_Peak(a) = max(mean(ccAltRes_excAc_forEstimate{a})) ; % variance peaks
    acAltInhRes_forEstimate_Peak(a) = max(mean(ccAltRes_inhAc_forEstimate{a})) ;
    
    gainExc(a) = max(mean(gAlt_Exc{a})) ; % gain factors from g
    gainInh(a) = max(mean(gAlt_Inh{a})) ;
    gainProduct(a) = gainExc(a)*gainInh(a) ;
end
gainExc_squared = gainExc.^2 ;
gainInh_squared = gainInh.^2 ;

temp = corrcoef(gainProduct,ccAltRes_forEstimate_Peak) ;
covEstimateQuality_linCoef = temp(1,2) ; % ability to estimate cov from gain product

temp = corrcoef(gainExc_squared,acAltExcRes_forEstimate_Peak) ;
varExcEstimateQuality_linCoef = temp(1,2) ; % ability to estimate var from gain

temp = corrcoef(gainInh_squared,acAltInhRes_forEstimate_Peak) ;
varInhEstimateQuality_linCoef = temp(1,2) ; % ability to estimate var from gain

fitCoefs = polyfit(gainProduct,ccAltRes_forEstimate_Peak,1) ;
CommonNoiseVar_estimateFromCov = fitCoefs(1) ; % best fit common noise variance from covariance estimate

fitCoefs = polyfit(gainExc_squared,acAltExcRes_forEstimate_Peak,1) ;
CommonNoiseVar_estimateFromVarExc = fitCoefs(1) ; % best fit common noise variance from variance estimate
%ExcNoiseVar_estimateFromVarExc = fitCoefs(2) ; % best fit of independant exc noise variance from variance estimate

fitCoefs = polyfit(gainInh_squared,acAltInhRes_forEstimate_Peak,1) ;
CommonNoiseVar_estimateFromVarInh = fitCoefs(1) ; % best fit common noise variance from variance estimate
%InhNoiseVar_estimateFromVarInh = fitCoefs(2) ; % best fit of independant exc noise variance from variance estimate

CommonNoiseVar_estimate = [CommonNoiseVar_estimateFromCov, CommonNoiseVar_estimateFromVarExc, CommonNoiseVar_estimateFromVarInh] ;

% CommonNoiseVar_estimate = (covEstimateQuality_linCoef*CommonNoiseVar_estimateFromCov + varExcEstimateQuality_linCoef*CommonNoiseVar_estimateFromVarExc + varInhEstimateQuality_linCoef*CommonNoiseVar_estimateFromVarInh)/...
%     (covEstimateQuality_linCoef + varExcEstimateQuality_linCoef + varInhEstimateQuality_linCoef) ;
% 
% corrCoefAltRes_estimate = gainProduct*CommonNoiseVar_estimate./...
%     sqrt((gainExc_squared.*CommonNoiseVar_estimate + ExcNoiseVar_estimateFromVarExc).*(gainInh_squared.*CommonNoiseVar_estimate + InhNoiseVar_estimateFromVarInh)) ;

% to save memmory clear unused variables
clearvars -except Monitor gAlt_Exc_ON_Mean gExc_ON_Mean gAlt_Inh_ON_Mean...
    gInh_ON_Mean gAlt_Exc_OFF_Mean gAlt_Inh_OFF_Mean gExc_OFF_Mean gInh_OFF_Mean...
    time_cc ccAltResMean ccResMean ccAltResShuffMean ccRes_SEM ccAltRes_SEM ccAltResShuff_SEM...
    time_ccON cc_ON ccAlt_ON time_ccOFF cc_OFF ccAlt_OFF ccAltRes_ON_Mean...
    ccAltResShuff_ON_Mean ccAltResShuff_ON_Mean ccAltRes_OFF_Mean...
    ccAltResShuff_OFF_Mean Ba gAlt_Exc_ON_max_Mean gAlt_Inh_ON_max_Mean...
    gExc_ON_max_Mean gInh_ON_max_Mean gAlt_Exc_OFF_max_Mean gAlt_Inh_OFF_max_Mean...
    gExc_OFF_max_Mean gInh_OFF_max_Mean NumBars SI A Input fp SpatialStimParams...
    frameRate ccRes_ON_Mean ccRes_OFF_Mean gAlt_Exc gAlt_Inh gAlt_Exc_ON gExc_ON gAlt_Inh_ON...
    gInh_ON gAlt_Exc_OFF gAlt_Inh_OFF gExc_OFF gInh_OFF gAlt_Exc_Res gAlt_Inh_Res...
    id2 id ccAltInd_ON_Mean ccAltInd_OFF_Mean time ccRes ccAltRes_Peak ccRes_Peak...
    gAlt_Exc_max_Mean gAlt_Inh_max_Mean gAlt_Exc_max_std gAlt_Inh_max_std...
    residualOption gAlt_Exc_Mean gAlt_Inh_Mean gainProduct gainExc_squared gainInh_squared...
    ccAltRes_forEstimate_Peak acAltExcRes_forEstimate_Peak acAltInhRes_forEstimate_Peak...
    covEstimateQuality_linCoef varExcEstimateQuality_linCoef varInhEstimateQuality_linCoef...
    CommonNoiseVar_estimate ccAltRes_Peak peakVector
    


% cell attached analysis
epochs = str2num(Input(A).(id2)) ;
for a = 1:length(epochs) ; % for each spike epoch
    [dataCA(a,:), error] = ITCReadEpoch(epochs(a), 0, fp) ;    % get data
        
    [SI_CA(a), error] = ITCGetSamplingInterval(epochs(a), fp); % get sampling interval
    SI_CA(a) = SI_CA(a) * 1e-6; % Sampling interval in sec
end

if Input(A).ITC18flag == 1 ;
    SI_CA = SI_CA*1.25 ;
end


time_CA = [1:length(dataCA)]*SI_CA(1) ;

for a = 1:length(epochs) ;
    StrucString = ['params_epoch',num2str(epochs(a))] ; 
    Struct = SpatialStimParams.(StrucString) ;

    BarAngles_CA(a,:) = Struct.BarAngle ;
    BarWidth_CA(a) = Struct.BarWidth ;
    BarSpeed_CA(a) = Struct.BarSpeed ;
    
    prePnts_CA(a) = floor(Struct.spatial_prepts/(frameRate*SI(1))) ;
    postPnts_CA(a) = floor(Struct.spatial_postpts/(frameRate*SI(1))) ;
    stimPnts_CA(a) = floor(Struct.spatial_stimpts/(frameRate*SI(1))) ;
    
    [Monitor_CA(a,:), error] = ITCReadEpoch(epochs(a), 1, fp) ;    % get data
    
    NumBars_CA(a) = length(BarAngles_CA(a,:)) ; % number of bars shown per trial
    BarPnts_CA(a) = floor(stimPnts_CA(a)/NumBars_CA(a)) ; % number of points the bar is presented + interbar points
    OnPnts_CA(a) = floor((BarWidth_CA(a)/BarSpeed_CA(a))/(frameRate*SI_CA(1))) ; % number of points before end of bar appears triggering off response
end

SpikePnts = SpikeDetection(dataCA,10,1/SI_CA(1)) ;
SpikeTrain = zeros(size(dataCA)) ;

for a=1:length(SpikePnts) ;
    SpikeTrain(a,SpikePnts{a}) = 1 ;
end

[Ba_CA,i] = sort(BarAngles_CA,2) ; % sort each row of bar angle in accending order
for  a = 1:NumBars_CA(1) ; % for each bar
    for b= 1:size(dataCA,1) ; % on each trial

        SpikeTrainBar{a}(b,:) = SpikeTrain(b,prePnts_CA(b)+BarPnts_CA(b)*(i(b,a)-1):prePnts_CA(b)+BarPnts_CA(b)*(i(b,a))) ;
        SpikeTrain_ON{a}(b,:) = SpikeTrain(b,prePnts_CA(b)+BarPnts_CA(b)*(i(b,a)-1):prePnts_CA(b)+BarPnts_CA(b)*(i(b,a)-1)+OnPnts_CA(b)) ;
        SpikeTrain_OFF{a}(b,:) = SpikeTrain(b,prePnts_CA(b)+BarPnts_CA(b)*(i(b,a)-1)+OnPnts_CA(b):prePnts_CA(b)+BarPnts_CA(b)*(i(b,a))) ;
        
        timeCheck{a}(b,:) = time_CA(1,prePnts_CA(b)+BarPnts_CA(b)*(i(b,a)-1):prePnts_CA(b)+BarPnts_CA(b)*(i(b,a))) ;
        timeCheck_ON{a}(b,:) = time_CA(1,prePnts_CA(b)+BarPnts_CA(b)*(i(b,a)-1):prePnts_CA(b)+BarPnts_CA(b)*(i(b,a)-1)+OnPnts_CA(b)) ;
        timeCheck_OFF{a}(b,:) = time_CA(1,prePnts_CA(b)+BarPnts_CA(b)*(i(b,a)-1)+OnPnts_CA(b):prePnts_CA(b)+BarPnts_CA(b)*(i(b,a))) ;
        
        CAexample{a}(b,:) = dataCA(b,prePnts_CA(b)+BarPnts_CA(b)*(i(b,a)-1):prePnts_CA(b)+BarPnts_CA(b)*(i(b,a))) ; % cell attached example
    
        FR{a}(b) = (sum(SpikeTrainBar{a}(b,:))/length(SpikeTrainBar{a}(b,:)))/SI_CA(b) ;
        FR_ON{a}(b) = (sum(SpikeTrain_ON{a}(b,:))/length(SpikeTrain_ON{a}(b,:)))/SI_CA(b) ;
        FR_OFF{a}(b) = (sum(SpikeTrain_OFF{a}(b,:))/length(SpikeTrain_OFF{a}(b,:)))/SI_CA(b) ;
    end
    
    MeanFR(a) = mean(FR{a}) ;
    MeanFR_ON(a) = mean(FR_ON{a}) ;
    MeanFR_OFF(a) = mean(FR_OFF{a}) ;
    
    StdFR(a) = std(FR{a}) ;
    StdFR_ON(a) = std(FR_ON{a}) ;
    StdFR_OFF(a) = std(FR_OFF{a}) ;
end


% figure
% figure
% subplot(3,1,1)
% plot([1:length(Monitor)],Monitor)
% xlabel('sample pnts')
% ylabel('Monitor reading')
% title('conductances')
% 
% subplot(3,1,2)
% plot([1:length(Monitor_CA)],Monitor_CA)
% xlabel('sample pnts')
% ylabel('Monitor reading')
% title('cell attached')

% subplot(3,1,3)
% imagesc(droppedFrame) ;


% figure
% 
% set(gcf,'position',[101 25 1439 1064])
% for a=1:NumBars ;
%     subplot(8,2,a*2-1)
%     plot([1:length(gAlt_Exc_ON_Mean{a})]*SI(1),gAlt_Exc_ON_Mean{a})
%     hold on
%     plot([1:length(gExc_ON_Mean{a})]*SI(1),gExc_ON_Mean{a},'b--')
%     plot([1:length(gAlt_Inh_ON_Mean{a})]*SI(1),gAlt_Inh_ON_Mean{a},'r')
%     plot([1:length(gInh_ON_Mean{a})]*SI(1),gInh_ON_Mean{a},'r--')
% %     xlabel('time')
% %     ylabel('conductance (nS)')
%     if a==1 ;
%         title('ON alt and sh conductances')
%     end
%     
%     subplot(8,2,a*2)
%     plot([1:length(gAlt_Exc_OFF_Mean{a})]*SI(1),gAlt_Exc_OFF_Mean{a})
%     hold on
%     plot([1:length(gExc_OFF_Mean{a})]*SI(1),gExc_OFF_Mean{a},'b--')
%     plot([1:length(gAlt_Inh_OFF_Mean{a})]*SI(1),gAlt_Inh_OFF_Mean{a},'r')
%     plot([1:length(gInh_OFF_Mean{a})]*SI(1),gInh_OFF_Mean{a},'r--')
% %     xlabel('time')
% %     ylabel('conductance (nS)')
%     if a==1 ;
%         title('OFF alt and sh conductances')
%     end
%     
% end
% legend('alt exc', 'sh exc', 'alt inh', 'sh inh')
% 
% figure
% set(gcf,'position',[101 25 1439 1064])
% for a=1:NumBars ;
%     subplot(8,4,a*4-3)
%     plot([1:length(gExc_ON_Mean{a})]*SI(1),gExc_ON{a},'b')
%     hold on
%     plot([1:length(gAlt_Exc_ON_Mean{a})]*SI(1),gAlt_Exc_ON{a},'k')
%    
%     subplot(8,4,a*4-2)
%     plot([1:length(gInh_ON_Mean{a})]*SI(1),gInh_ON{a},'r')
%     hold on
%     plot([1:length(gAlt_Inh_ON_Mean{a})]*SI(1),gAlt_Inh_ON{a},'k')
% %     xlabel('time')
% %     ylabel('conductance (nS)')
%     if a==1 ;
%         title('ON alt and sh conductances')
%     end
%     
%     subplot(8,4,a*4-1)
%     plot([1:length(gExc_OFF_Mean{a})]*SI(1),gExc_OFF{a},'b')    
%     hold on
%     plot([1:length(gAlt_Exc_OFF_Mean{a})]*SI(1),gAlt_Exc_OFF{a},'k')
%     
%     subplot(8,4,a*4)    
%     plot([1:length(gInh_OFF_Mean{a})]*SI(1),gInh_OFF{a},'r')
%     hold on
%     plot([1:length(gAlt_Inh_OFF_Mean{a})]*SI(1),gAlt_Inh_OFF{a},'k')    
% %     xlabel('time')
% %     ylabel('conductance (nS)')
%     if a==1 ;
%         title('OFF alt and sh conductances')
%     end
%     
% end
% legend('alt exc', 'sh exc', 'alt inh', 'sh inh')
% 
% 
% figure
% set(gcf,'position',[101 25 1439 1064])
% for a=1:NumBars ;
%     subplot(8,2,a*2-1)
%     plot(time_ccON,cc_ON(a,:),'color',[1 0 1])
%     hold on
%     plot(time_ccON,ccAlt_ON(a,:),'k')
%     plot(time_ccON,ccAltInd_ON_Mean(a,:),'k--')
%     if a==1 ;
%         title('ON alt and sh cc of mean g')
%     end
%     
%     subplot(8,2,a*2)
%     plot(time_ccOFF,cc_OFF(a,:),'color',[1 0 1])
%     hold on
%     plot(time_ccOFF,ccAlt_OFF(a,:),'k')
%     plot(time_ccOFF,ccAltInd_OFF_Mean(a,:),'k--')
%     if a==1 ;
%         title('OFF alt and sh cc of mean g')
%     end
% end
% legend('sh cc','alt cc','alt ind cc')
% 
% figure
% set(gcf,'position',[101 25 1439 1064])
% for a=1:NumBars ;
%     subplot(8,2,a*2-1)
%     plot(time_ccON,ccAltRes_ON_Mean(a,:),'k')
%     hold on
%     plot(time_ccON,ccAltResShuff_ON_Mean(a,:),'g')
%     if a==1 ;
%         title('ON alt and shuffled cc of res g')
%     end
%     
%     subplot(8,2,a*2)
%     plot(time_ccOFF,ccAltRes_OFF_Mean(a,:),'k')
%     hold on
%     plot(time_ccOFF,ccAltResShuff_OFF_Mean(a,:),'g')
%     if a==1 ;
%         title('OFF alt and shuffled cc of res g')
%     end
% 
% end
% legend('alt','shuff')
% 
% figure
% set(gcf,'position',[101 25 1439 1064])
% subplot(4,2,1)
% errorbar(Ba_CA(1,:),MeanFR_ON,-StdFR_ON,StdFR_ON,'c*-')
% xlabel('bar angle')
% ylabel('mean firing rate')
% title('ON response')
% set(gca,'xlim',[0 360]) 
% 
% subplot(4,2,2)
% errorbar(Ba_CA(1,:),MeanFR_OFF,-StdFR_OFF,StdFR_OFF,'c*-')
% xlabel('bar angle')
% ylabel('mean firing rate')
% title('OFF response')
% set(gca,'xlim',[0 360]) 
% 
% subplot(4,2,3)
% plot(Ba(1,:),gAlt_Exc_ON_max_Mean,'b*-')
% hold on
% plot(Ba(1,:),gAlt_Inh_ON_max_Mean,'r*-')
% plot(Ba(1,:),gExc_ON_max_Mean,'b*--')
% plot(Ba(1,:),gInh_ON_max_Mean,'r*--')
% xlabel('bar angle')
% ylabel('peak conductance mean (nS)')
% legend('alt exc','alt inh','sh exc', 'sh inh')
% set(gca,'xlim',[0 360])
% 
% subplot(4,2,4)
% plot(Ba(1,:),gAlt_Exc_OFF_max_Mean,'b*-')
% hold on
% plot(Ba(1,:),gAlt_Inh_OFF_max_Mean,'r*-')
% plot(Ba(1,:),gExc_OFF_max_Mean,'b*--')
% plot(Ba(1,:),gInh_OFF_max_Mean,'r*--')
% xlabel('bar angle')
% ylabel('peak conductance mean (nS)')
% legend('alt exc','alt inh','sh exc', 'sh inh')
% set(gca,'xlim',[0 360])
% 
% subplot(4,2,5)
% plot(Ba(1,:),ccAlt_ON(:,time_ccON==0),'k*-')
% hold on
% plot(Ba(1,:),cc_ON(:,time_ccON==0),'*-','color',[1 0 1])
% plot(Ba(1,:),ccAltRes_ON_Mean(:,time_ccON==0),'k*--')
% plot(Ba(1,:),ccRes_ON_Mean(:,time_ccON==0),'*--','color',[1 0 1])
% xlabel('bar angle')
% ylabel('cc')
% legend('alt mean cc','sh mean cc','alt res cc','sh res cc')
% set(gca,'xlim',[0 360])
% 
% subplot(4,2,6)
% plot(Ba(1,:),ccAlt_OFF(:,time_ccOFF==0),'k*-')
% hold on
% plot(Ba(1,:),cc_OFF(:,time_ccOFF==0),'*-','color',[1 0 1])
% plot(Ba(1,:),ccAltRes_OFF_Mean(:,time_ccOFF==0),'k*--')
% plot(Ba(1,:),ccRes_OFF_Mean(:,time_ccOFF==0),'*--','color',[1 0 1])
% xlabel('bar angle')
% ylabel('cc')
% legend('alt mean cc','sh mean cc','alt res cc','sh res cc') 
% set(gca,'xlim',[0 360])
%     
% subplot(4,2,7)
% plot(MeanFR_ON,gExc_ON_max_Mean/max(gExc_ON_max_Mean),'b*')
% hold on
% plot(MeanFR_ON,gInh_ON_max_Mean/max(gInh_ON_max_Mean),'r*')
% plot(MeanFR_ON,cc_ON(:,time_ccON==0),'k*')
% plot(MeanFR_ON,gAlt_Exc_ON_max_Mean/max(gAlt_Exc_ON_max_Mean),'bo')
% plot(MeanFR_ON,gAlt_Inh_ON_max_Mean/max(gAlt_Inh_ON_max_Mean),'ro')
% plot(MeanFR_ON,ccAlt_ON(:,time_ccON==0),'ko')
% xlabel('firing rate')
% ylabel('exc, inh, cc norm')
% legend('sh exc','sh inh','sh cc','alt exc','alt inh','alt cc') 
% 
% subplot(4,2,8)
% plot(MeanFR_OFF,gExc_OFF_max_Mean/max(gExc_OFF_max_Mean),'b*')
% hold on
% plot(MeanFR_OFF,gInh_OFF_max_Mean/max(gInh_OFF_max_Mean),'r*')
% plot(MeanFR_OFF,cc_OFF(:,time_ccOFF==0),'k*')
% plot(MeanFR_OFF,gAlt_Exc_OFF_max_Mean/max(gAlt_Exc_OFF_max_Mean),'bo')
% plot(MeanFR_OFF,gAlt_Inh_OFF_max_Mean/max(gAlt_Inh_OFF_max_Mean),'ro')
% plot(MeanFR_OFF,ccAlt_OFF(:,time_ccOFF==0),'ko')
% 
% xlabel('firing rate')
% ylabel('exc, inh, cc norm')
% legend('sh exc','sh inh','sh cc','alt exc','alt inh','alt cc')  

% figure
% set(gcf,'position',[101 25 1439 1064])
% for a=1:8
%     subplot(8,1,a)
%     errorbar(time_cc,ccAltResMean(a,:),-ccAltRes_SEM(a,:),ccAltRes_SEM(a,:),'k')
%     hold on
%     errorbar(time_cc,ccAltResShuffMean(a,:),-ccAltResShuff_SEM(a,:),ccAltResShuff_SEM(a,:),'g')
%     errorbar(time_cc,ccResMean(a,:),-ccRes_SEM(a,:),ccRes_SEM(a,:),'r')
% end

% figure
% subplot(2,1,1)
% errorbar(Ba(1,:),gAlt_Exc_max_Mean,-gAlt_Exc_max_std,gAlt_Exc_max_std)
% hold on
% errorbar(Ba(1,:),gAlt_Inh_max_Mean,-gAlt_Inh_max_std,gAlt_Inh_max_std,'r')
% xlabel('bar direction')
% ylabel('G peak')
% 
% subplot(2,1,2)
% errorbar(Ba(1,:),ccAltRes_Peak,'k')
% hold on
% plot(Ba(1,:),ccRes_Peak,'g')
% xlabel('bar direction')
% ylabel('cc peak')
% 
% figure
% plot(MeanFR/max(MeanFR),ccAltRes_Peak,'k*')
% hold on 
% plot(MeanFR/max(MeanFR),ccRes_Peak,'g*')
% 
% figure
% errorbar(Ba_CA(1,:),MeanFR, StdFR, StdFR) 

figure  % this figure add 10/15/10
%set(gcf,'position',[101 25 1439 1064])
for a=1:8
    subplot(1,8,a)
    plot(time_cc,ccAltResMean(a,:),'k')
    hold on
    plot(time_cc,ccResMean(a,:),'g')
    set(gca,'xlim',[-.2,.2],'ylim',[-.5,1])
end

figure
plot(gainProduct,ccAltRes_forEstimate_Peak,'*')

figure
plot(gainExc_squared,acAltExcRes_forEstimate_Peak,'*')
hold on
plot(gainInh_squared,acAltInhRes_forEstimate_Peak,'r*')


% for igor 
% 
%mean cross correlations for each bar direction
for a=1:8 ;
    identifier = ['ccAltResMean',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = ccAltResMean(a,:) ;
    
%     identifier = ['ccAltResShuffMean',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = ccAltResShuffMean(a,:) ;
    
%     identifier = ['ccResMean',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = ccResMean(a,:) ;
    
end

identifier = ['time_cc','cell',num2str(A)] ;
ForIgor.(identifier) = time_cc ;

% idividual alternating conductances, residuals and res cross correlations
concat_gAlt_Exc = cell2mat(gAlt_Exc) ; % concatinate bar directions
concat_gAlt_Inh = cell2mat(gAlt_Inh) ;

concat_gAlt_ExcRes = cell2mat(gAlt_Exc_Res) ;
concat_gAlt_InhRes = cell2mat(gAlt_Inh_Res) ;

maxExc = max(max(concat_gAlt_Exc)) ;
maxInh = max(max(concat_gAlt_Inh)) ;


for a=1:size(concat_gAlt_Exc,1) ; 
    
%     % normalized for g clamp
%     identifier = ['ExcG1t',num2str(a),'c',num2str(A)] ;
%     ForIgor.(identifier) = concat_gAlt_Exc(a,:)/maxExc ;
% 
%     identifier = ['InhG1t',num2str(a),'c',num2str(A)] ;
%     ForIgor.(identifier) = concat_gAlt_Inh(a,:)/maxInh ;

    
    % unormalized
    identifier = ['exGaltExc',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = concat_gAlt_Exc(a,:) ;

    identifier = ['exGaltInh',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = concat_gAlt_Inh(a,:) ;

    identifier = ['exGaltExcRes',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = concat_gAlt_ExcRes(a,:) ;
    
    identifier = ['exGaltInhRes',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = concat_gAlt_InhRes(a,:) ;    
end

identifier = ['time','cell',num2str(A)] ;
ForIgor.(identifier) = [1:length(concat_gAlt_Exc)]*SI(1) ;
 
% % mean cross corr coef for each bar angle
% for a=1:8 ;
%     identifier = ['corrCoefAltResMean',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = ccAltResMean(a,time_cc==0) ;
%     
%     identifier = ['corrCoefAltResShuffMean',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = ccAltResShuffMean(a,time_cc==0) ;   
% end
% 
% % cell attached mean firing rate for polar plots
% 
% identifier = ['meanFR',num2str(A)] ;
% ForIgor.(identifier) = [MeanFR,MeanFR(1)] ;
% 
% identifier = ['barAngles',num2str(A)] ;
% ForIgor.(identifier) = [Ba_CA(1,:),Ba_CA(1,1)] ;


% % cell attached mean firing rate and peak of cc
% identifier = ['meanFRCA',num2str(A)] ;
% ForIgor.(identifier) = MeanFR ;
% 
% identifier = ['meanFRCAnorm',num2str(A)] ;
% ForIgor.(identifier) = MeanFR/max(MeanFR) ;
% 
% identifier = ['barAnglesLin',num2str(A)] ;
% ForIgor.(identifier) = Ba_CA(1,:) ;

% % cell attached examples
% 
% for b=1:NumBars ;
%     stCA = nans(size(SpikeTrainBar{b})) ;
%     stCA(SpikeTrainBar{b}==1)=1 ;
% 
%     for a=1:size(stCA,1) ; % for each trial
%         identifier = ['SpikeTrainCA',num2str(b),'trial',num2str(a),'cell',num2str(A)] ;
%         ForIgor.(identifier) = stCA(a,:)*a ;
%     end
% end
% 
% identifier = ['timeCAexample',num2str(A)] ;
% ForIgor.(identifier) = timeCheck{1}(1,:) ;

% 
% identifier = ['ccAltResPeak',num2str(A)] ;
% ForIgor.(identifier) = ccAltRes_Peak ;
% 
% identifier = ['ccResPeak',num2str(A)] ;
% ForIgor.(identifier) = ccRes_Peak ;
% 
% % G exc/inh tuning curves
% 
% identifier = ['meanGExcpeak',num2str(A)] ;
% ForIgor.(identifier) = gAlt_Exc_max_Mean ;
% 
% identifier = ['meanGInhpeak',num2str(A)] ;
% ForIgor.(identifier) = gAlt_Inh_max_Mean ;
% 
% identifier = ['meanGExcstd',num2str(A)] ;
% ForIgor.(identifier) = gAlt_Exc_max_std ;
% 
% identifier = ['meanGInhstd',num2str(A)] ;
% ForIgor.(identifier) = gAlt_Inh_max_std ;
% 
% identifier = ['barAnglesLin',num2str(A)] ;
% ForIgor.(identifier) = Ba(1,:) ;

% %mean cross correlations for each bar direction using optional residual calc
% if residualOption==1 ;
%     for a=1:8 ;    
%         identifier = ['ccAltResMeanResOpt',num2str(a),'cell',num2str(A)] ;
%         ForIgor.(identifier) = ccAltResMean(a,:) ;
% 
%         identifier = ['ccResMeanResOpt',num2str(a),'cell',num2str(A)] ;
%         ForIgor.(identifier) = ccResMean(a,:) ;
% 
%     end
% end

% estimating noise correlation based on tuning curves and simple model
identifier = ['covExcInh','cell',num2str(A)] ;
ForIgor.(identifier) = ccAltRes_forEstimate_Peak ;

identifier = ['varExc','cell',num2str(A)] ;
ForIgor.(identifier) = acAltExcRes_forEstimate_Peak ;

identifier = ['varInh','cell',num2str(A)] ;
ForIgor.(identifier) = acAltInhRes_forEstimate_Peak ;

identifier = ['gainProduct','cell',num2str(A)] ;
ForIgor.(identifier) = gainProduct ;

identifier = ['gainExcSquared','cell',num2str(A)] ;
ForIgor.(identifier) = gainExc_squared ;

identifier = ['gainInhSquared','cell',num2str(A)] ;
ForIgor.(identifier) = gainInh_squared ;

identifier = ['CovEstimateLinCorr','cell',num2str(A)] ;
ForIgor.(identifier) = covEstimateQuality_linCoef ;

identifier = ['varExcEstimateLinCorr','cell',num2str(A)] ;
ForIgor.(identifier) = varExcEstimateQuality_linCoef ;

identifier = ['varInhEstimateLinCorr','cell',num2str(A)] ;
ForIgor.(identifier) = varInhEstimateQuality_linCoef ;

identifier = ['NcEstimates','cell',num2str(A)] ;
ForIgor.(identifier) = CommonNoiseVar_estimate ;

% mean conductances
concat_gAlt_Exc_Mean = cell2mat(gAlt_Exc_Mean) ; % concatinate bar directions
concat_gAlt_Inh_Mean = cell2mat(gAlt_Inh_Mean) ;

identifier = ['meanGaltExc','cell',num2str(A)] ;
ForIgor.(identifier) = concat_gAlt_Exc_Mean ;

identifier = ['meanGaltInh','cell',num2str(A)] ;
ForIgor.(identifier) = concat_gAlt_Inh_Mean ;
 
% tuning curve of peak cc aligned by peak to a 1:8 vector peak at 4
identifier = ['ccPeakVector','cell',num2str(A)] ;
ForIgor.(identifier) = peakVector ;

%ForIgor.nada = 1 ;


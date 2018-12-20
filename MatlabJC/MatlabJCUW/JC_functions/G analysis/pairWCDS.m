function ForIgor = pairWCDS(Input,Parameters,id,A) ; 
% this function will analyze wc data from a pair of ds cells in responses to a moving
% bar

% 10/13/10 edditted to calculate local residual

for input = 1:length(id) ; % for each input

% get data
[fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/mouse/',Input(A{input}).cellname]);

epochs = str2num(Input(A{input}).(id{input})) ;
numTrials = length(epochs) ;

for a = 1:numTrials ; 
    [data1(a,:), error] = ITCReadEpoch(epochs(a), 0, fp) ;    % get data
    [data2(a,:), error] = ITCReadEpoch(epochs(a), 1, fp) ;    % get data
    
    [SI(a), error] = ITCGetSamplingInterval(epochs(a), fp); % get sampling interval
    SI(a) = SI(a) * 1e-6; % Sampling interval in sec
end

if Input(A{input}).ITC18flag == 1 ;
    SI = SI*1.25 ;
end
time = [1:length(data1)]*SI(1) ;

% low pass filter to remove electrical crap
data1_lpf = lowPassFilter(data1, 1/SI(1), 5000) ; %(signal,samplerate,cutoff frequ (hz))
data2_lpf = lowPassFilter(data2, 1/SI(1), 5000) ; 

% high pass filter to remove slow drift 
% greg S. sent this to me to help implement a butterworth and avoid ringing
F=1 ; % filter cuttoff
Wn = F*SI(1); %normalized frequency cutoff
[z, p, k] = butter(1,Wn,'high'); %
[sos,g]=zp2sos(z,p,k); 
myfilt=dfilt.df2sos(sos,g);

data1_hpf = filter(myfilt,data1_lpf')'; % filter implementation
data2_hpf = filter(myfilt,data2_lpf')'; 

% subtract off means
for a=1:numTrials ;
    data1_hpf_meansub(a,:) = data1_hpf(a,:) - mean(data1_hpf(a,:)) ;
    data2_hpf_meansub(a,:) = data2_hpf(a,:) - mean(data2_hpf(a,:)) ;
end

% change to conductances 
if strcmp(id{input},'ExcDS') ;
    G1 = data1_hpf_meansub/-62 ;
    G2 = data2_hpf_meansub/-62 ;
elseif strcmp(id{input},'InhDS') ;
    G1 = data1_hpf_meansub/62 ;
    G2 = data2_hpf_meansub/62 ;
elseif strcmp(id{input},'InhExcDS') ;
    G1 = data1_hpf_meansub/62 ;
    G2 = data2_hpf_meansub/-62 ;
elseif strcmp(id{input},'ExcInhDS') ;
    G1 = data1_hpf_meansub/-62 ;
    G2 = data2_hpf_meansub/62 ;
end

% offset
G1 = G1-min(G1(:)) ;
G2 = G2-min(G2(:)) ;

% get stimulus parameters
SpatialStimParams = hdf5load(['~/Data/mouse/',Input(A{input}).cellname,'_spatial.h5']) ; % load spatial stim params
frameRate = 60 ;


for a = 1:numTrials ;
    StrucString = ['params_epoch',num2str(epochs(a))] ; 
    Struct = SpatialStimParams.(StrucString) ;

    BarAngles(a,:) = Struct.BarAngle ;
    BarWidth(a) = Struct.BarWidth ;
    BarSpeed(a) = Struct.BarSpeed ;
    
    prePnts(a) = floor(Struct.spatial_prepts/(frameRate*SI(1))) ;
    postPnts(a) = floor(Struct.spatial_postpts/(frameRate*SI(1))) ;
    groupPnts(a) = floor((length(data1)-prePnts(a)-postPnts(a))/length(BarAngles(a,:))) ;
    
    [Monitor(a,:), error] = ITCReadEpoch(epochs(a), 5, fp) ;    % get data
    
    NumBars(a) = length(BarAngles(a,:)) ; % number of bars shown per trial
    BarPnts(a) = floor((length(data1)-prePnts(a)-postPnts(a))/NumBars(a)) ; % number of points the bar is presented + interbar points
    OnPnts(a) = floor((BarWidth(a)/BarSpeed(a))/(frameRate*SI(1))) ; % number of points before end of bar appears triggering off response
end


% get G and G stats for each bar direction

[Ba,i] = sort(BarAngles,2) ; % sort each row of bar angle in accending order
for  a = 1:NumBars(1) ; % for each bar
    for b= 1:numTrials ; % on each trial get spike number

        GBar1{a}(b,:) = G1(b,prePnts(b)+BarPnts(b)*(i(b,a)-1):prePnts(b)+BarPnts(b)*(i(b,a))) ;
        GBar1_ON{a}(b,:) = G1(b,prePnts(b)+BarPnts(b)*(i(b,a)-1):prePnts(b)+BarPnts(b)*(i(b,a)-1)+OnPnts(b)) ;
        GBar1_OFF{a}(b,:) = G1(b,prePnts(b)+BarPnts(b)*(i(b,a)-1)+OnPnts(b):prePnts(b)+BarPnts(b)*(i(b,a))) ;
        
        GBar2{a}(b,:) = G2(b,prePnts(b)+BarPnts(b)*(i(b,a)-1):prePnts(b)+BarPnts(b)*(i(b,a))) ;
        GBar2_ON{a}(b,:) = G2(b,prePnts(b)+BarPnts(b)*(i(b,a)-1):prePnts(b)+BarPnts(b)*(i(b,a)-1)+OnPnts(b)) ;
        GBar2_OFF{a}(b,:) = G2(b,prePnts(b)+BarPnts(b)*(i(b,a)-1)+OnPnts(b):prePnts(b)+BarPnts(b)*(i(b,a))) ;    
        
        timeCheck{a}(b,:) = time(1,prePnts(b)+BarPnts(b)*(i(b,a)-1):prePnts(b)+BarPnts(b)*(i(b,a))) ;
        timeCheck_ON{a}(b,:) = time(1,prePnts(b)+BarPnts(b)*(i(b,a)-1):prePnts(b)+BarPnts(b)*(i(b,a)-1)+OnPnts(b)) ;
        timeCheck_OFF{a}(b,:) = time(1,prePnts(b)+BarPnts(b)*(i(b,a)-1)+OnPnts(b):prePnts(b)+BarPnts(b)*(i(b,a))) ;    
        
    end
    
    GBar1_mean(a,:) = mean(GBar1{a}) ;
    GBar1_mean_ON(a,:) = mean(GBar1_ON{a}) ;
    GBar1_mean_OFF(a,:) = mean(GBar1_OFF{a}) ;
    
    GBar2_mean(a,:) = mean(GBar2{a}) ;
    GBar2_mean_ON(a,:) = mean(GBar2_ON{a}) ;
    GBar2_mean_OFF(a,:) = mean(GBar2_OFF{a}) ;
    
    % residuals
    residualOption = 1 ;
    if residualOption == 0 ; % get residuals standard calculation
        for b = 1:numTrials ;
            GBar1_res{a}(b,:) = GBar1{a}(b,:) -  GBar1_mean(a,:) ;
            GBar2_res{a}(b,:) = GBar2{a}(b,:) -  GBar2_mean(a,:) ;
        end

    elseif residualOption == 1 ; % get residual nearest neighbor calculation (added option 8/4/10)
        for b= 2:numTrials-1 ;
            GBar1_res{a}(b-1,:) = GBar1{a}(b,:) -  (GBar1{a}(b-1,:)+GBar1{a}(b+1,:))/2 ;
            GBar2_res{a}(b-1,:) = GBar2{a}(b,:) -  (GBar1{a}(b-1,:)+GBar1{a}(b+1,:))/2 ;
        end
    end
    
    
    
    for b= 1:numTrials-residualOption*2 ; % residuals and cc
        cc{a}(b,:) = xcov(GBar1_res{a}(b,:),GBar2_res{a}(b,:)) ;
        cc_coef{a}(b,:) = xcov(GBar1_res{a}(b,:),GBar2_res{a}(b,:),'coef') ;
        GBar1_res_ac(b) = max(xcov(GBar1_res{a}(b,:))) ; % auttocorr peak
        GBar2_res_ac(b) = max(xcov(GBar2_res{a}(b,:))) ;
    end

    cc_mean(a,:) = mean(cc{a}) ;
    cc_coef_mean(a,:) = mean(cc_coef{a}) ; 
    
    cc_peak(input,a) = CCpeakFinder(cc_mean(a,:)) ; % peak of cc 
    cc_coef_peak(input,a) = CCpeakFinder(cc_coef_mean(a,:)) ; % peak of cc coef 
    
    GBar1_res_ac_peak(input,a) = mean(GBar1_res_ac) ;
    GBar2_res_ac_peak(input,a) = mean(GBar2_res_ac) ;
        
    % peak g
    GBar1_mean_peak(a)= max(GBar1_mean(a,:)) ;
    GBar1_mean_peak_ON(a)= max(GBar1_mean_ON(a,:)) ;
    GBar1_mean_peak_OFF(a)= max(GBar1_mean_OFF(a,:)) ;   
    
    GBar2_mean_peak(a)= max(GBar2_mean(a,:)) ;
    GBar2_mean_peak_ON(a)= max(GBar2_mean_ON(a,:)) ;
    GBar2_mean_peak_OFF(a)= max(GBar2_mean_OFF(a,:)) ;
    
end % bar angle loop

timeCC=SI(1)*([1:length(cc_mean)]-(length(cc_mean)+1)/2) ;


% % check stimulus precision
% MonitorSwing=nans(numTrials,3*ceil(max(time))) ;
% for a=1:numTrials ;
%     Monitor_norm(a,:) = Monitor(a,:)-min(Monitor(a,:)) ;
%     Monitor_norm(a,:) = Monitor_norm(a,:)/max(Monitor_norm(a,:)) ;
%     Monitor_norm(a,:) = smooth(Monitor_norm(a,:),5) ;
%     temp = time(find(diff(Monitor_norm(a,:)>.5)==1)) ;
%     MonitorSwingNum(a) = length(temp) ;
%     MonitorSwing(a,1:length(temp)) = temp ;
% end
% 
% [r,c] = find(diff(MonitorSwing)>0.010) ;

%figures

figure
plot(time,Monitor)
hold on
if ~isempty(r)
    plot(time,Monitor(r+1,:),'r--')
    title(num2str(unique(r')+1))
end

% figure % raw data divided by bar cell 1
% for a=1:NumBars(1) ;
%     subplot(NumBars(1)/2,2,a)
%     plot(timeCheck{1}(1,:),GBar1{a})
% end
% 
% figure % raw data divided by bar cell 2
% for a=1:NumBars(1) ;
%     subplot(NumBars(1)/2,2,a)
%     plot(timeCheck{1}(1,:),GBar2{a})
% end
% 
% for a=1:NumBars(1) ; % residula data divided into figure by bar and trial by pannel 
%     figure
%     for b=1:numTrials ;
%         subplot(ceil(numTrials/4),4,b)
%         plot(timeCheck{1}(1,:),GBar1_res{a}(b,:)-mean(GBar1_res{a}(b,:)),'b')
%         hold on
%         plot(timeCheck{1}(1,:),GBar2_res{a}(b,:)-mean(GBar2_res{a}(b,:)),'r')
%     end
%     
% end


figure % peak of mean current
plot(Ba,GBar1_mean_peak,'b*-')
hold on
plot(Ba,GBar2_mean_peak,'r*-')
xlabel('bar angle')
ylabel('peak of mean g')
legend('cell 1','cell 2')

figure
for a=1:NumBars(1) ;
    subplot(NumBars(1)/2,2,a)
    plot(timeCC,cc_coef_mean(a,:))
end

clearvars -except cc_coef_peak cc_peak GBar1_res_ac_peak GBar2_res_ac_peak ...
    Input A id

end % input loop

% calculate bounds for simulataneous measure
Et1 = [] ;
It1 = [] ;
Et2 = [] ;
It2 = [] ;
Et1Et2 = [] ;
It1It2 = [] ;
Et1It2 = [] ;
Et2It1 = [] ;

for a=1:length(id) ;
    if strcmp(id{a},'ExcInhDS') ;
        Et1It2 = cc_peak(a,:) ;
        Et1 = [Et1;GBar1_res_ac_peak(a,:)] ;
        It2 = [It2;GBar2_res_ac_peak(a,:)] ;
        
    elseif strcmp(id{a},'InhExcDS') ;
        Et2It1 = cc_peak(a,:) ;
        Et2 = [Et2;GBar2_res_ac_peak(a,:)] ;
        It1 = [It1;GBar1_res_ac_peak(a,:)] ;
        
    elseif strcmp(id{a},'ExcDS') ;
        Et1Et2 = cc_peak(a,:) ;
        Et1 = [Et1;GBar1_res_ac_peak(a,:)] ;
        Et2 = [Et2;GBar2_res_ac_peak(a,:)] ;
        
    elseif strcmp(id{a},'InhDS') ;
        It1It2 = cc_peak(a,:) ;
        It1 = [It1;GBar1_res_ac_peak(a,:)] ;
        It2 = [It2;GBar2_res_ac_peak(a,:)] ;
        
    end
end

Et1 = mean(Et1) ; % average over all estimates
It1 = mean(It1) ;
Et2 = mean(Et2) ;
It2 = mean(It2) ;

if ~isempty(Et1Et2) && ~isempty(It1It2)
    
    for a=1:length(Et1) ;
        WeakUpperBound_Et1It1(a) = (min(Et1Et2(a),It1It2(a)) + min(Et1(a)-Et1Et2(a),It1(a)-It1It2(a)))/sqrt(Et1(a)*It1(a)) ;
        WeakUpperBound_Et2It2(a) = (min(Et1Et2(a),It1It2(a)) + min(Et2(a)-Et1Et2(a),It2(a)-It1It2(a)))/sqrt(Et2(a)*It2(a)) ;
    end
end

if ~isempty(Et1Et2) && ~isempty(It1It2) && ~isempty([Et1It2,Et2It1])

    EIc = mean([Et1It2;Et2It1]) ;
         
    for a=1:length(Et1) ;
        UpperBound_Et1It1(a) = (EIc(a) + min(Et1(a)-Et1Et2(a),It1(a)-It1It2(a)))/sqrt(Et1(a)*It1(a)) ;
        UpperBound_Et2It2(a) = (EIc(a) + min(Et2(a)-Et1Et2(a),It2(a)-It1It2(a)))/sqrt(Et2(a)*It2(a)) ;
        
        LowerBound(a) = EIc(a)/sqrt(Et1(a)*It1(a)) ;
    end
    
    Estimate = EIc./sqrt(Et1Et2(a)*It1It2(a)) ;
end



figure
plot(UpperBound_Et1It1,'b*-')
hold on
plot(UpperBound_Et2It2,'r*-')
plot(LowerBound,'k*-')






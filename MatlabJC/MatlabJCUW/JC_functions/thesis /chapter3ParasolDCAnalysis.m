function ForIgor = chapter3ParasolDCAnalysis(Input,A)

% thesis chapter 3 parasol dynamic clamp experiments

% see cells 070908Bc1,092408Ac2,"c3,"Bc1,"Bc2,"Bc3, injected with 1st set of g 
% see cells 102108Bc1,"c2,"c3,"c4,"c5,"c6,102208Bc1,"c2 injected with 2nd

% example figure
preTime = 1 ; %sec from light stimulus
postTime = .1 ; % sec

[fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/primate/',Input(A).cellname]);
 
epochsControl = str2num(Input(A).dcControl) ; 
epochsMinusInh = str2num(Input(A).dcMinusInh) ;

for a = 1:length(epochsControl) ; % for each g epoch
    [dataControl(a,:),error] = ITCReadEpoch(epochsControl(a), 0, fp) ;    % get data
    [excgControl(a,:), inhgControl(a,:), error] = ITCReadEpochStmGClamp(epochsControl(a), 0, fp) ;
end

for a = 1:length(epochsMinusInh) ; % for each g epoch
    [dataMinusInh(a,:),error] = ITCReadEpoch(epochsMinusInh(a), 0, fp) ;
    [excgMinusInh(a,:), inhgMinusInh(a,:), error] = ITCReadEpochStmGClamp(epochsMinusInh(a), 0, fp) ;
end

[SI, error] = ITCGetSamplingInterval(epochsControl(1), fp); % get sampling interval
SI = SI*10^-6 ;

prePnts = round(preTime/SI) ;
postPnts = round(postTime/SI) ;

time = [1:length(dataControl)]*SI ;

% spike detection
spikePntsControl = SpikeDetection_WC(dataControl,-20,1/SI) ;
spikePntsMinusInh = SpikeDetection_WC(dataMinusInh,-20, 1/SI) ;

spikeTrainControl = zeros(size(dataControl)) ;
spikeTrainMinusInh = zeros(size(dataMinusInh)) ;
for a=1:length(epochsControl) ;
    spikeTrainControl(a,spikePntsControl{a}) = 1 ;
end
for a=1:length(epochsMinusInh) ;
    spikeTrainMinusInh(a,spikePntsMinusInh{a}) = 1 ;
end

% number of spikes
spikeNumControl = sum(spikeTrainControl(:,prePnts:end-postPnts),2) ;
spikeNumMinusInh = sum(spikeTrainMinusInh(:,prePnts:end-postPnts),2) ;

for a=1:5 ;
    spikeNumMeanControl(a) = mean(spikeNumControl(a:5:end)) ; 
    spikeNumMeanMinusInh(a) = mean(spikeNumMinusInh(a:5:end)) ;
    
    spikeNumSEMControl(a) = std(spikeNumControl(a:5:end))/sqrt(length(spikeNumControl(a:5:end))) ; 
    spikeNumSEMMinusInh(a) = std(spikeNumMinusInh(a:5:end))/sqrt(length(spikeNumMinusInh(a:5:end))) ; 
end
    
% voltage histograms
voltageBins = [-100:1:20] ;
VhistControl = hist(dataControl(:),voltageBins) ;
VhistMinusInh = hist(dataMinusInh(:),voltageBins) ;

VpdfControl = VhistControl/sum(VhistControl) ;
VpdfMinusInh = VhistMinusInh/sum(VhistMinusInh) ;

% % spike threshold 2nd deriv
% SpikeThreshold1Control = SpikeThresholdFinder(dataControl, spikePntsControl, SI, 1) ;
% SpikeThreshold1MinusInh = SpikeThresholdFinder(dataMinusInh, spikePntsMinusInh, SI, 1) ;
% 
% SpikeThreshold1HistControl = hist(cell2mat(SpikeThreshold1Control),voltageBins) ;
% SpikeThreshold1HistMinusInh = hist(cell2mat(SpikeThreshold1MinusInh),voltageBins) ;
% 
% SpikeThreshold1PdfControl = SpikeThreshold1HistControl/sum(SpikeThreshold1HistControl) ;
% SpikeThreshold1PdfMinusInh = SpikeThreshold1HistMinusInh/sum(SpikeThreshold1HistMinusInh) ;

% spike threshold peak voltage
SpikeThreshold2Control = SpikeThresholdFinder(dataControl, spikePntsControl, SI, 2) ;
SpikeThreshold2MinusInh = SpikeThresholdFinder(dataMinusInh, spikePntsMinusInh, SI, 2) ;

SpikeThreshold2HistControl = hist(SpikeThreshold2Control,voltageBins) ;
SpikeThreshold2HistMinusInh = hist(SpikeThreshold2MinusInh,voltageBins) ;

SpikeThreshold2PdfControl = SpikeThreshold2HistControl/sum(SpikeThreshold2HistControl) ;
SpikeThreshold2PdfMinusInh = SpikeThreshold2HistMinusInh/sum(SpikeThreshold2HistMinusInh) ;

SpikeThreshold2meanControl = mean(SpikeThreshold2Control) ;
SpikeThreshold2meanMinusInh = mean(SpikeThreshold2MinusInh) ;

% sta
load('~/data_analysis/LightStimParasolDC') ; % light stimulus
staTime = .3 ; %sec
staLength = staTime/SI ; % pnts
staControl = staFinder(repmat(LightStimParasolDC(:,prePnts:end-postPnts),length(epochsControl)/5,1),spikeTrainControl(:,prePnts:end-postPnts),staLength) ;
staMinusInh = staFinder(repmat(LightStimParasolDC(:,prePnts:end-postPnts),length(epochsMinusInh)/5,1),spikeTrainMinusInh(:,prePnts:end-postPnts),staLength) ;

staTime = [-staLength+1:0]*SI ;


% figures

figure
plot(time,dataControl(1,:),'k')
hold on
plot(time,dataMinusInh(1,:),'r')

figure
errorbar(spikeNumMeanControl,spikeNumMeanMinusInh,spikeNumSEMMinusInh)
hold on
plot(spikeNumMeanControl,spikeNumMeanControl,'k-')
xlabel('spike num control')
ylabel('spike num minus inh')

figure
plot(voltageBins,VpdfControl)
hold on
plot(voltageBins,VpdfMinusInh,'r')
xlabel('voltage bins')
ylabel('probability')

% figure
% plot(voltageBins,SpikeThreshold1PdfControl)
% hold on
% plot(voltageBins,SpikeThreshold1PdfMinusInh,'r')
% xlabel('spike threshold bins')
% ylabel('probability')

figure
plot(voltageBins,SpikeThreshold2PdfControl,'--')
hold on
plot(voltageBins,SpikeThreshold2PdfMinusInh,'r--')
xlabel('spike threshold bins')
ylabel('probability')

figure
plot(staTime,staControl)
hold on
plot(staTime,staMinusInh,'r')
xlabel('time')
ylabel('sta')

% For Igor

% example dynamic clamp voltage
identifier = ['VdataControlExample',num2str(A)] ;
ForIgor.(identifier) = dataControl(1,:) ; 

identifier = ['VdataMinusInhExample',num2str(A)] ;
ForIgor.(identifier) = dataMinusInh(1,:) ; 

identifier = ['VdataExampleTime',num2str(A)] ;
ForIgor.(identifier) = time ; 

% example spike trains
spikeTrainControl(1,(spikeTrainControl(1,:)==0)) = nan ;
spikeTrainMinusInh(1,(spikeTrainMinusInh(1,:)==0)) = nan ;

identifier = ['strainControlExample',num2str(A)] ;
ForIgor.(identifier) = spikeTrainControl(1,:) ; 

identifier = ['strainMinusInhExample',num2str(A)] ;
ForIgor.(identifier) = spikeTrainMinusInh(1,:)*2 ; 

% exmple conductances
identifier = ['gexcControlExample',num2str(A)] ;
ForIgor.(identifier) = excgControl(1,:) ; 

identifier = ['ginhControlExample',num2str(A)] ;
ForIgor.(identifier) = inhgControl(1,:) ;

identifier = ['gexcMinusInhExample',num2str(A)] ;
ForIgor.(identifier) = excgMinusInh(1,:) ; 

identifier = ['ginhMinusInhExample',num2str(A)] ;
ForIgor.(identifier) = inhgMinusInh(1,:) ;

% voltage histograms
identifier = ['VhistBins',num2str(A)] ;
ForIgor.(identifier) = voltageBins ; 

identifier = ['VpdfControl',num2str(A)] ;
ForIgor.(identifier) = VpdfControl ; 

identifier = ['VpdfMinusInh',num2str(A)] ;
ForIgor.(identifier) = VpdfMinusInh ; 

% spike threshold
identifier = ['sThresh2pdfControl',num2str(A)] ;
ForIgor.(identifier) = SpikeThreshold2PdfControl ; 

identifier = ['sThresh2pdfMinusInh',num2str(A)] ;
ForIgor.(identifier) = SpikeThreshold2PdfMinusInh ; 

identifier = ['sThreshMeanControl',num2str(A)] ;
ForIgor.(identifier) = SpikeThreshold2meanControl ; 

identifier = ['sThreshMeanMinusInh',num2str(A)] ;
ForIgor.(identifier) = SpikeThreshold2meanMinusInh ; 

% sta
identifier = ['staTime',num2str(A)] ;
ForIgor.(identifier) = staTime ;

identifier = ['staControl',num2str(A)] ;
ForIgor.(identifier) = staControl ;

identifier = ['staMinusInh',num2str(A)] ;
ForIgor.(identifier) = staMinusInh ;





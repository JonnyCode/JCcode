function ForIgor = CAnoiseAnalyzer(Input,Parameters,id,A) ;

% this will analyze dynamic clamp data from midget cells injected with a 
gaussianStd = 10.^[-3:.1:-1] ;

epochs = str2num(Input(A).(id)) ;

[fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/primate/',Input(A).cellname]);

for a = 1:length(epochs) ;
    [Data(a,:), error] = ITCReadEpoch(epochs(a), 0, fp) ; 
    [stm(a,:), error] = ITCReadEpochStm(epochs(a), 0, fp);  % get the light stimulus
end

[SI, error] = ITCGetSamplingInterval(epochs(1), fp); % get sampling interval    
SI = SI * 1e-6; % Sampling interval in sec
if Input(A).ITC18flag == 1 ;
    SI = SI*1.25 ;
end
SI = .0001 ;
time = [SI:SI:SI*length(Data)] ;

[prePnts, error] = ITCGetStmPrePts(epochs(1), 0, 0, fp) ; % points collected beyond which data is not worth analyzing
[postPnts, error] = ITCGetStmTailPts(epochs(1), 0, 0, fp) ;

SpikePnts = SpikeDetection(Data,25,10000) ;
SpikeTrain = zeros(size(Data)) ;
for a=1:length(epochs) ;
    SpikeTime{a} = time(SpikePnts{a}) ;
    SpikeTrain(a,SpikePnts{a}) = 1 ;
end

[meanPSTHsnr,sumPSTHvar,psth] = PsthVar(SpikeTrain,gaussianStd ,1/SI) ; 
minIntTime = interp1(meanPSTHsnr,gaussianStd,1) ;
psthi = find(gaussianStd>.005,1,'first') ;

[powerX,meanSpikeSpectrum,resSpikeSpectrum,meanSpikeSpectrum_smth,resSpikeSpectrum_smth,snrSpikeSpectrum_smth] = snrSpikeSpectrum(SpikeTrain,1/SI,.001) ;

firingRate = sum(SpikeTrain(:,prePnts:end-postPnts),2)/(time(end-postPnts)-time(prePnts)) ;


figure
rasterPlot(SpikeTime)
xlabel('time (sec)')
ylabel('trial')

figure
plot(gaussianStd,meanPSTHsnr)
xlabel('gaussian std')
ylabel('mean psth snr')


figure
plot(time,mean(psth{psthi}))
xlabel('time (sec)')
ylabel('firing rate (hz)')

ForIgor.time = time ;
ForIgor.MeanPSTH = mean(psth{psthi}) ;
ForIgor.firingRate = firingRate ;
ForIgor.gaussianStd = gaussianStd ;
ForIgor.meanPSTHsnr = meanPSTHsnr ;
ForIgor.minIntTime = minIntTime ;
for a=1:length(epochs); iden = ['spikeTrain',num2str(a)], ForIgor.(iden) = SpikeTrain(a,:) ; end






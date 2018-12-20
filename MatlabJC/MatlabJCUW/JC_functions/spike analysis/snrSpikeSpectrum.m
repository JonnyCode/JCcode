function [powerX,meanSpikeSpectrum,resSpikeSpectrum,meanSpikeSpectrum_smth,resSpikeSpectrum_smth,snrSpikeSpectrum_smth,VarSumResSpikeSpectrum] = snrSpikeSpectrum(spikeTrain,samplerate,spikeTimeResolution) 

% JC '09
% spikeTrain (1s and 0s), samplerate (hz), spikeTimeResolution - time over which you can resolve a single spike
% (not counting refractory period, sec) 
% JC edit '10 to add unbiased estimate SNR (van Hateren and Snippe, 2001)

NumTrains = size(spikeTrain,1) ;
spikeTrainIntStep = samplerate*spikeTimeResolution ;

% create spike trains at requested resolution
for a=1:NumTrains ;
    newSpikeTrain(a,:) = DecimateWave(spikeTrain(a,:),spikeTrainIntStep)*spikeTrainIntStep ; 
end

% mean spike train
meanSpikeTrain = mean(newSpikeTrain) ;

% get residual spike trains and noise spectrum
resSpikeTrain = newSpikeTrain - repmat(meanSpikeTrain,size(newSpikeTrain,1),1) ;

for a=1:NumTrains ;
    [powerX,resSpikeSpectrum_ind(a,:)] = PowerSpectrumFinder(resSpikeTrain(a,:),1/spikeTimeResolution) ; %power spectrum of individual residuals
end
resSpikeSpectrum_biased = mean(resSpikeSpectrum_ind) ;
resSpikeSpectrum = resSpikeSpectrum_biased*NumTrains/(NumTrains-1) ; %correcting bias

% get signal spectrum
[powerX,meanSpikeSpectrum] = PowerSpectrumFinder(meanSpikeTrain,1/spikeTimeResolution) ; %power spectrum
meanSpikeSpectrum = meanSpikeSpectrum - resSpikeSpectrum_biased/(NumTrains-1) ; % correcting bias

% get varaince of noise amplitude between trials
SumResSpikeSpectrumInd = sum(resSpikeSpectrum_ind,2) ; % noise var on each trial
VarSumResSpikeSpectrum = var(SumResSpikeSpectrumInd) ; % variance of noise var between trials

% smooth power spectrums
meanSpikeSpectrum_smth = SmoothPowerSpectrum(meanSpikeSpectrum(2:end), powerX(2:end), 2, 1) ; % smoothed power spectrum 
resSpikeSpectrum_smth = SmoothPowerSpectrum(resSpikeSpectrum(2:end), powerX(2:end), 2, 1) ;

% calculate unbiased snr from smoothed spectrums
snrSpikeSpectrum_smth = meanSpikeSpectrum_smth.PowerSpec./resSpikeSpectrum_smth.PowerSpec ; % snr of smoothed power spectrum



% figure
% subplot(4,1,1)
% plot(powerX,meanSpikeSpectrum,'k*-')
% hold on
% plot(meanSpikeSpectrum_smth.Freq,meanSpikeSpectrum_smth.PowerSpec,'b-') ;
% h=gca;,set(h,'xscale','log','yscale','log') ;
% subplot(4,1,2)
% plot(powerX,resSpikeSpectrum,'k*-')
% hold on
% plot(resSpikeSpectrum_smth.Freq,resSpikeSpectrum_smth.PowerSpec,'r-') ;
% h=gca;,set(h,'xscale','log','yscale','log') ;
% subplot(4,1,3)
% plot(meanSpikeSpectrum_smth.Freq,meanSpikeSpectrum_smth.PowerSpec,'b-') ;
% hold on
% plot(resSpikeSpectrum_smth.Freq,resSpikeSpectrum_smth.PowerSpec,'r-') ;
% h=gca;,set(h,'xscale','log','yscale','log') ;
% subplot(4,1,4)
% plot(resSpikeSpectrum_smth.Freq,snrSpikeSpectrum_smth,'k-') ;
% hold on
% plot(resSpikeSpectrum_smth.Freq,1,'r-') ;
% h=gca;,set(h,'xscale','log') ;
end


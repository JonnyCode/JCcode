
% continuous case
pnts = 50000 ;
samplerate = 10000 ;

signal = normrnd(1,1,1,pnts) ;
gaussianStd = 10.^[-3:.1:-1] ;
frequs = [samplerate,1000,100,10] ;
trialNumber = 10 ;

for a =1:length(frequs)
    for b=1:trialNumber
        Noise = normrnd(0,.25,1,pnts) ;
        Noise = lowPassFilter(Noise,samplerate,frequs(a)) ;
        input{a}(b,:) = signal + Noise ;
    end
    [meanPSTHsnr(a,:),sumPSTHvar(a,:),psth{a}] = PsthVar(input{a},gaussianStd,samplerate) ;
    [powerX,meanSpikeSpectrum,resSpikeSpectrum,meanSpikeSpectrum_smth,resSpikeSpectrum_smth,snrSpikeSpectrum_smth(a,:)] = snrSpikeSpectrum(input{a},samplerate,1/samplerate) ; 
end

time = [1/samplerate:1/samplerate:(1/samplerate)*pnts] ;
for a=1:length(gaussianStd) ;
    gaussian(a,:) = exp(-((time-(time(end)/2)).^2)/(2*gaussianStd(a)^2)) ;
    [powerspec_xvalues(a,:), mean_powerspec(a,:)] = PowerSpectrumFinder(gaussian(a,:)/sum(gaussian(a,:)),samplerate);
end

for a=1:length(gaussianStd) ;
    HalfPower = max(mean_powerspec(a,:))*.5 ;
    i=find(mean_powerspec(a,:)<HalfPower,1,'first');
    GaussianRolloff(a) = powerspec_xvalues(a,i) ;
end
    

figure
plot(time,gaussian)

figure
plot(powerspec_xvalues(1,:), mean_powerspec)
h=gca ;
set(h,'xScale','log')

figure
plot(gaussianStd,GaussianRolloff,'*-')

figure
plot(gaussianStd,sumPSTHvar,'*')
h=gca ;
set(h,'xScale','log')

figure
plot(gaussianStd,meanPSTHsnr,'*')
h=gca ;
set(h,'xScale','log')

figure
plot(GaussianRolloff, meanPSTHsnr,'*-')

figure
plot(meanSpikeSpectrum_smth.Freq,snrSpikeSpectrum_smth)



% after conversion to spike train
%filter = 
NFgain = .2 ;
NFthresh = 0 ;
NFmax = 1 ;
spikeThresh = .7 ;
ref = .002 ;


for a = 1:length(frequs) ;
    spikeTrain{a} = zeros(trialNumber,length(time)) ;
    for b=1:trialNumber
        if exist('filter','var')
            inputLinear{a}(b,:) = conv(input{a}(b,:),filter) ; % linear filter
        else
            inputLinear{a}(b,:)= input{a}(b,:) ;
        end
        %nonlinear filter
        inputNLpre = inputLinear{a}(b,:)*NFgain - NFgain*NFthresh ; % gain
        inputNLpre(inputLinear{a}(b,:)<NFthresh) = 0 ; %threshold   
        inputNLpre(inputNLpre>NFmax) = NFmax ; % saturate

        inputNL{a}(b,:) = inputNLpre(1,1:length(time)) ; % cut off convolution extra
        
        %spike generation
        forspikes = inputNL{a}(b,:) ; 
        i = find(inputNL{a}(b,:)>spikeThresh,1,'first');
        if ~isempty(i)
        while ~isempty(i) ;
            spikeTrain{a}(b,i)=1;
            i = find(forspikes>spikeThresh,1,'first') ;
            forspikes(i:i+samplerate*ref)=0 ;
        end
        spiketimes{a}{b}=time(spikeTrain{a}(b,:)==1) ;
        end 
    end
    
    [meanPSTHsnr(a,:),sumPSTHvar(a,:),psth{a}] = PsthVar(spikeTrain{a},gaussianStd,samplerate) ;
    [powerX,meanSpikeSpectrum,resSpikeSpectrum,meanSpikeSpectrum_smth(a),resSpikeSpectrum_smth(a),snrSpikeSpectrum_smth(a,:)] = snrSpikeSpectrum(spikeTrain{a},samplerate,1/samplerate) ; 
    
    figure
    rasterPlot(spiketimes{a})
end

figure
plot(gaussianStd,sumPSTHvar,'*')
h=gca ;
set(h,'xScale','log')

figure
plot(gaussianStd,meanPSTHsnr,'*')
h=gca ;
set(h,'xScale','log')

       
figure
plot(meanSpikeSpectrum_smth(1).Freq,snrSpikeSpectrum_smth)
set(gca, 'xscale','log') 

figure
for a=1:length(frequs);
plot(meanSpikeSpectrum_smth(a).Freq,meanSpikeSpectrum_smth(a).PowerSpec)
hold on
plot(resSpikeSpectrum_smth(a).Freq,resSpikeSpectrum_smth(a).PowerSpec,'r')
end
set(gca, 'xscale','log') 

        
        
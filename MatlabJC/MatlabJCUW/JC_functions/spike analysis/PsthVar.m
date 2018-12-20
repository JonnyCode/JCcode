function [meanPSTHsnr,sumPSTHvar,psth] = PsthVar(spikeTrain,gaussianStd,samplerate) 

% gaussianStd (sec), samplerate (hz),

time = [1/samplerate:1/samplerate:(1/samplerate)*length(spikeTrain)] ; 
convfactor = floor(length(time)/2) ;

for a=1:length(gaussianStd) ; % for each new psth window gaussian
    gaussian = exp(-((time-(time(end)/2)).^2)/(2*gaussianStd(a)^2)) ; % gaussian function exp(-((t-u)^2/2s^2)) from wikipedia "gauassian functions"
    
    for b = 1:size(spikeTrain,1)  % for each set of spike trains       
        tempConv = ifft(fft([spikeTrain(b,:),zeros(1,length(time))]).*fft([gaussian,zeros(1,length(time))])); % convolution w/ zero padding to avoid inapproprate wrap around      
        psthTemp= tempConv(convfactor:end-convfactor)*samplerate/sum(gaussian) ; % conv zero padds at begining and end (thus convfactor) and normalize to firing rate
        %psthTemp=tempConv(convfactor:end-convfactor)/sum(tempConv(convfactor:end-convfactor)) ; % this line normalizes by sum of psth   
        psth{a}(b,:) =psthTemp(1:length(time)) ; % corrects length if spikeTrain is odd length
    end
    
    sumPSTHvar(a) = sum(var(psth{a},[],1)) ;
    
    meanPSTH(a,:) = mean(psth{a},1) ;
    stdPSTH(a,:) = std(psth{a},[],1) ;
    
    meanPSTHsnr(a)=mean(meanPSTH(a,meanPSTH(a,:)>0)./stdPSTH(a,meanPSTH(a,:)>0)) ;
    
end

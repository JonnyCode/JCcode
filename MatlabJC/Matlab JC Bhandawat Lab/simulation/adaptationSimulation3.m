% simulation to test equation designed to maximinze response entropy and
% minimize response noise (i.e. maximize mutual information).  Different
% from adapatationSimulation2 because I will be testing nongaussian
% distributions (and Field and Rieke only).

% questions: what is the optimal NL to increase discriminability between
% these signal and noise distributions.  Should the nonlinearity be based
% only on the signal or on the total distribution containing signal and
% noise?  Does dividing by the summ of the nonlinearities help?

% choose response distribution type
InputRspRangeMin = -200 ;
InputRspRangeDelta = .1 ;
InputRspRangeMax = 200 ;

% number of inputs to average over
InputNum = 20 ;
InputSampleNumber = 10000000 ; % number of input samples per dist (num of samples in output dist will be this /InputNum)

% parameters
InputRspRange = [InputRspRangeMin:InputRspRangeDelta:InputRspRangeMax] ;
InputRspDistMean = [39, 45] ;
InputRspDistStd = [20, 20] ;
FracSignal = .5 ; 

% sampled data
InputSample = normrnd(InputRspDistMean(1),InputRspDistStd(1),2,InputSampleNumber) ;
InputSample_sig = normrnd(InputRspDistMean(2),InputRspDistStd(2),1,floor(InputSampleNumber*FracSignal)) ;
InputSample(2,1:length(InputSample_sig)) = InputSample_sig ;
randindex = randperm(InputSampleNumber) ;
InputSample(2,:) = InputSample(2,randindex) ;

% Gaussian response distributions
InputRspDist = nan(length(InputRspDistMean),length(InputRspRange)) ;
for a=1:length(InputRspDistMean) ;
    InputRspDist(a,:) = hist(InputSample(a,:),InputRspRange) ; 
end

InputRspDist_sig = hist(InputSample_sig,InputRspRange) ;

% nonlinearity
NL1 = InputRspRange.*(InputRspDist(2,:)./sum(InputRspDist,1)) ; % NL based on entire distributions
NL2 = InputRspRange.*(InputRspDist_sig./(InputRspDist(1,:)+InputRspDist_sig)) ; % NL based only on signal hiden in noise

NotNani = ~isnan(NL1) ;
NL1 = NL1(NotNani) ;
NL1x = InputRspRange(NotNani) ;

NotNani = ~isnan(NL2) ;
NL2 = NL2(NotNani) ;
NL2x = InputRspRange(NotNani) ;

% output distributions

% nonlinearity
for a = 1:length(InputRspDistMean) ;
    InputSample_NL1(a,:) = interp1(NL1x,NL1,InputSample(a,:),'nearest') ; % apply nonlinearity to sample data
    InputSample_NL2(a,:) = interp1(NL2x,NL2,InputSample(a,:),'nearest') ;
end

% summation
for a = 1:length(InputRspDistMean) ;
    c = 0;
    for b=1:InputNum:InputSampleNumber ;
        c = c+1 ;
        OutputSamples_linear(a,c) = mean(InputSample(a,b:b+InputNum-1)) ;
        OutputSamples_NL1(a,c) = mean(InputSample_NL1(a,b:b+InputNum-1)) ;
        OutputSamples_NL2(a,c) = mean(InputSample_NL2(a,b:b+InputNum-1)) ;
    end
end

OutputRspDist_linear = nan(length(InputRspDistMean),length(InputRspRange)) ;
OutputRspDist_NL1 = nan(length(InputRspDistMean),length(InputRspRange)) ;
for a=1:length(InputRspDistMean) ;
    OutputRspDist_linear(a,:) = hist(OutputSamples_linear(a,:),InputRspRange) ;
    OutputRspDist_NL1(a,:) = hist(OutputSamples_NL1(a,:),InputRspRange) ;
    OutputRspDist_NL2(a,:) = hist(OutputSamples_NL2(a,:),InputRspRange) ;
end

% probability transform
for a = 1:length(InputRspDistMean) ;
    InputSample_pt1(a,:) = interp1(InputRspRange,InputRspDist(1,:)./sum(InputRspDist(1,:)),InputSample(a,:),'linear') ;
    InputSample_pt2(a,:) = interp1(InputRspRange,InputRspDist(2,:)./sum(InputRspDist(2,:)),InputSample(a,:),'linear') ;
end
for a = 1:length(InputRspDistMean) ;
    c = 0;
    for b=1:InputNum:InputSampleNumber ;
        c = c+1 ;
        OutputSamples_pt(a,c) = prod(InputSample_pt2(a,b:b+InputNum-1))/(prod(InputSample_pt1(a,b:b+InputNum-1))+prod(InputSample_pt2(a,b:b+InputNum-1))) ;
    end
end
for a=1:length(InputRspDistMean) ;
    OutputRspDist_pt(a,:) = hist(OutputSamples_pt(a,:),InputRspRange) ;
end


% discriminnability
DiscrimInput = sum(max(InputRspDist,[],1))/sum(sum(InputRspDist,1)) 
DiscrimOutput_linear = sum(max(OutputRspDist_linear,[],1))/sum(sum(OutputRspDist_linear,1)) 
DiscrimOutput_NL1 = sum(max(OutputRspDist_NL1,[],1))/sum(sum(OutputRspDist_NL1,1)) 
DiscrimOutput_NL2 = sum(max(OutputRspDist_NL2,[],1))/sum(sum(OutputRspDist_NL2,1))
DiscrimOuput_pt = sum(max(OutputRspDist_pt,[],1))/sum(sum(OutputRspDist_pt,1))

% figures
figure
subplot(5,1,1)
plot(InputRspRange,InputRspDist)

subplot(5,1,2)
plot(InputRspRange,OutputRspDist_linear)

subplot(5,1,3)
plot(InputRspRange,OutputRspDist_NL1)

subplot(5,1,4)
plot(InputRspRange,OutputRspDist_NL2)

subplot(5,1,5)
plot(InputRspRange,OutputRspDist_pt)

figure
plotyy(InputRspRange,InputRspDist,NL1x,NL1)



figure
plot(InputSample',InputSample_NL1','.')
hold on
plot(NL1x,NL1,'k-')
plot(NL1x,NL1x,'y-')


figure
subplot(2,2,4)
plot(InputRspRange,InputRspDist)
xlim([min(InputRspRange)-1,max(InputRspRange)+1])
subplot(2,2,1)
plot(OutputRspDist_linear,InputRspRange)
ylim([min(InputRspRange)-1,max(InputRspRange)+1])
subplot(2,2,2)
plot(InputRspRange,InputRspRange,'*-')
xlim([min(InputRspRange)-1,max(InputRspRange)+1])
ylim([min(InputRspRange)-1,max(InputRspRange)+1])
% simulation to test equation designed to maximinze response entropy and
% minimize response noise (i.e. maximize mutual information).  Different
% from adapatationSimulation1 in that I have abondoned attempt to deal send
% pdfs through nonlinearity and integration and work with simulated data
% instead.

% questions: does equation behave like field and rieke with discreet high
% noise distibutions and laughlin-barlow at continuous low noise
% distributions?  How does it behave in intermediates?  Does it improve
% linear discriminability?

% choose response distribution type
InputRspDistType = 'DiscreteNoisy' ;
InputRspRangeMin = -200 ;
InputRspRangeDelta = .1 ;
InputRspRangeMax = 200 ;

% number of inputs to average over
InputNum = 2000 ;
InputSampleNumber = 10000000 ; % number of input samples per dist (num of samples in output dist will be this /InputNum)

% parameters
InputRspRange = [InputRspRangeMin:InputRspRangeDelta:InputRspRangeMax] ;
if strcmp(InputRspDistType,'DiscreteNoisy') ; % discrete paramters (noisy distributions)
    InputRspDistAmp = [1,1] ;
    InputRspDistMean = [39, 39] ;
    InputRspDistStd = [20, 20] ;
elseif strcmp(InputRspDistType,'ContNoiseless') ; % pseudo continuous paramters (noiseless distributions)
    InputRspDistAmp =  exp(-((InputRspRange-InputRspRange(round(length(InputRspRange)/2))).^2)/(2*InputRspRange(round(length(InputRspRange)/4))^2))  ; 
    InputRspDistMean = [0:InputRspRangeDelta:100] ;
    InputRspDistStd = ones(1,length(InputRspRange))*10^-100 ;
elseif strcmp(InputRspDistType,'ContNoisyAdditive') ; % pseudo continuous paramters (noisy distributions addative)
    InputRspDistAmp =  exp(-((InputRspRange-InputRspRange(round(length(InputRspRange)/2))).^2)/(2*InputRspRange(round(length(InputRspRange)/4))^2))  ; 
    InputRspDistMean = [0:InputRspRangeDelta:100] ;
    InputRspDistStd = ones(1,length(InputRspRange))*50 ;
elseif strcmp(InputRspDistType,'ContNoisyMultiplicative') ; % pseudo continuous paramters (noisy distributions multiplicative)
    InputRspDistAmp =  exp(-((InputRspRange-InputRspRange(round(length(InputRspRange)/2))).^2)/(2*InputRspRange(round(length(InputRspRange)/4))^2))  ; 
    InputRspDistMean = [0:InputRspRangeDelta:100] ;
    InputRspDistStd = [1:length(InputRspRange)] ;
end

% sampled data
InputSample = nan(length(InputRspDistMean),InputSampleNumber) ;
for a=1:length(InputRspDistMean) ;
    InputSample(a,:) = normrnd(InputRspDistMean(a),InputRspDistStd(a),1,InputSampleNumber) ;
end

% Gaussian response distributions
InputRspDist = nan(length(InputRspDistMean),length(InputRspRange)) ;
for a=1:length(InputRspDistMean) ;
    InputRspDist(a,:) = hist(InputSample(a,:),InputRspRange) ; 
end

% intergrals
InputRspDist_indSum = nan(1,length(InputRspDistMean)) ;
for a=1:length(InputRspDistMean) ; %across responses within a single stimulus
    InputRspDist_indSum(a) = sum(InputRspDist(a,:)) ;
end

InputRspDist_sum = sum(InputRspDist,1) ; % across stimuli caused by a single response

% Oi (laughlin like - probability any response will come from given distribution - maximize entropy)
Oi = cumsum(InputRspDist_indSum)/sum(InputRspDist_indSum) ;

% Pi (rieke like - probability a particular response (r) will come from a given distribution (a) - maximize discriminability )
Pr = nan(length(InputRspDistMean),length(InputRspRange)) ;
for a = 1:length(InputRspDistMean) ;
    Pr(a,:) = InputRspDist(a,:)./InputRspDist_sum ;
end

% wieght Oi by Pi to maximize mutual info (I hope)
Or = nan(length(InputRspDistMean),length(InputRspRange)) ;
for a = 1:length(InputRspDistMean) ;
    Or(a,:) =  Oi(a) *Pr(a,:) ;
end

% nonlinearity
NL = sum(Or,1)*InputRspRangeMax ;

NL = InputRspRange.*(InputRspDist(2,:)./sum(InputRspDist,1)) ; % 2 dist try

NotNani = ~isnan(NL) ;
NL = NL(NotNani) ;
NLx = InputRspRange(NotNani) ;
% output distributions

% nonlinearity
for a = 1:length(InputRspDistMean) ;
    InputSample_NL(a,:) = interp1(NLx,NL,InputSample(a,:),'nearest') ; % apply nonlinearity to sample data
end

% summation
for a = 1:length(InputRspDistMean) ;
    c = 0;
    for b=1:InputNum:InputSampleNumber ;
        c = c+1 ;
        OutputSamples_linear(a,c) = mean(InputSample(a,b:b+InputNum-1)) ;
        OutputSamples_NL(a,c) = mean(InputSample_NL(a,b:b+InputNum-1)) ;
    end
end

OutputRspDist_linear = nan(length(InputRspDistMean),length(InputRspRange)) ;
OutputRspDist_NL = nan(length(InputRspDistMean),length(InputRspRange)) ;
for a=1:length(InputRspDistMean) ;
    OutputRspDist_linear(a,:) = hist(OutputSamples_linear(a,:),InputRspRange) ;
    OutputRspDist_NL(a,:) = hist(OutputSamples_NL(a,:),InputRspRange) ;
end


% discriminnability
DiscrimInput = sum(max(InputRspDist,[],1))/sum(sum(InputRspDist,1)) 
DiscrimOutput_linear = sum(max(OutputRspDist_linear,[],1))/sum(sum(OutputRspDist_linear,1)) 
DiscrimOutput_NL = sum(max(OutputRspDist_NL,[],1))/sum(sum(OutputRspDist_NL,1)) 


% figures
figure
subplot(3,1,1)
plot(InputRspRange,InputRspDist)

subplot(3,1,2)
plot(InputRspRange,OutputRspDist_linear)

subplot(3,1,3)
plot(InputRspRange,OutputRspDist_NL)

figure
plotyy(InputRspRange,InputRspDist,NLx,NL)



figure
plot(InputSample',InputSample_NL','.')
hold on
plot(InputRspRange,NL,'k-')
plot(InputSample(:),InputSample(:),'yo')


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

figure
subplot(2,2,4)
plot(InputRspRange,InputRspDist)
xlim([min(InputRspRange)-1,max(InputRspRange)+1])
subplot(2,2,1)
plot(OutputRspDist_NL,InputRspRange)
ylim([min(InputRspRange)-1,max(InputRspRange)+1])
subplot(2,2,2)
plot(InputRspRange,NL,'*-')
xlim([min(InputRspRange)-1,max(InputRspRange)+1])
ylim([min(InputRspRange)-1,max(InputRspRange)+1])



figure
subplot(3,1,1)
plot(InputRspDistMean,Oi,'*')

subplot(3,1,2)
plot(InputRspRange,Pr)

subplot(3,1,3)
plot(InputRspRange,Or)


figure
plot(InputRspRange,InputRspDist_sum,'b')
hold on
plot(OutputRspRange,OutputRspDist_sum,'r')

figure
plot(DiscrimInput)




% simple sim (delta Y seems to be much larger than delta X for this to work)
Xmax = 100 ;
deltaX = .001 ;
deltaY = 1 ;
Xrange = [0:deltaX:Xmax] ;
Yrange = [0:deltaY:Xmax] ;

x = exp(-((Xrange-50).^2)/(2*10^2)) ;
xx = x/sum(x) ;
y = cumsum(xx)*Xmax;

for a=1:length(Yrange)-1 ;
    i = find(y>=Yrange(a) & y<Yrange(a+1)) ;
    yy(a) = sum(xx(i)) ;
end

yy=yy/sum(yy) ;

figure
plot(Xrange,xx)
hold on
plot(Yrange(1:end-1),yy,'r')
% simulation to test equation designed to maximinze response entropy and
% minimize response noise (i.e. maximize mutual information)

% questions: does equation behave like field and rieke with discreet high
% noise distibutions and laughlin-barlow at continuous low noise
% distributions?  How does it behave in intermediates?  Does it improve
% linear discriminability?

% choose response distribution type
InputRspDistType = 'DiscreteNoisy' ;
InputRspRangeDelta = .01 ;
InputRspRangeMax = 100 ;

% parameters
InputRspRange = [0:InputRspRangeDelta:InputRspRangeMax] ;
if strcmp(InputRspDistType,'DiscreteNoisy') ; % discrete paramters (noisy distributions)
    InputRspDistAmp = [1,1,1] ;
    InputRspDistMean = [20, 40, 60] ;
    InputRspDistStd = [2, 2, 2] ;
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

% Gaussian response distributions
InputRspDist = nan(length(InputRspDistMean),length(InputRspRange)) ;
for a=1:length(InputRspDistMean) ;
    InputRspDist(a,:) = InputRspDistAmp(a)*exp(-((InputRspRange-InputRspDistMean(a)).^2)/(2*InputRspDistStd(a)^2)) ; 
end

% rectify response distributions
InputRspDist(InputRspDist>max(InputRspRange))= 0 ;
InputRspDist(InputRspDist<min(InputRspRange))= 0 ;

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

% % output distributions
% NL = NL + [10^-10:10^-10:length(NL)*10^-10] ; % this is a hack and should be removed once issue is resolved properly
% 
% OutputRspDist = nan(length(InputRspDistMean),length(InputRspRange)) ;
% for a=1:length(InputRspRange) ;
%     for b=1:length(InputRspDistMean) ;
%         if a==1 ;        
%             OutputRspDist(b,a) = InputRspDist(b,a)*InputRspRangeDelta/NL(a) ;
%             OutputRspRange(a) = NL(a) ;
%             
%         else ;
%             OutputRspDist(b,a) = ((InputRspDist(b,a-1)+InputRspDist(b,a))*InputRspRangeDelta/abs(NL(a-1)-NL(a)))-OutputRspDist(b,a-1) ;
%             OutputRspRange(a) = NL(a) ;
%         end
%     end
% end
    
%output distributions 
OutputRspDist = nan(length(InputRspDistMean),length(InputRspRange)) ;
OutputRspRange = NL ;
a=1 ;
b=1 ;
while a<=length(InputRspRange) ;
    if a==1 
        i = find(NL(a+1:end)~=NL(a),1,'first') ;
        if isempty(i) ;
            i = length(NL)-a+1 ;
        end
        if i>1 ;
            weightingFnct = [1,2*ones(1,i-1),1] ; 
        else
            weightingFnct = 1 ;
        end
        OutputRspDist(:,b) = (InputRspDist(:,a:a-1+i)*weightingFnct'*InputRspRangeDelta/NL(a)) ; 
        a = a + i  ;
        
    elseif a~=length(InputRspRange) ; 
        i = find(NL(a+1:end)~=NL(a),1,'first') ;
        if isempty(i) ;
            i = length(NL)-a+1 ;
        end
        weightingFnct = [1,2*ones(1,i-1),1] ; % NEEDS TO BE FIXED
        OutputRspDist(:,b) = (InputRspDist(:,a-1:a-1+i)*weightingFnct'*InputRspRangeDelta/(NL(a)-NL(a-1)))-OutputRspDist(:,b-1) ;
        a = a + i  ;
    end
    
    b=b+1 ;
end

OutputRspDist_sum = sum(OutputRspDist,1) ;

% discriminnability
for a=1:length(InputRspRange) ;
    DiscrimInput(a) = max(InputRspDist(:,a))/InputRspDist_sum(a) ;
end

for a=1:length(InputRspRange) ;
    DiscrimOutput(a) = max(OutputRspDist(:,a))/OutputRspDist_sum(a) ;
end


% figures
figure
plot([1:length(OutputRspDist)],OutputRspDist)


figure
subplot(2,2,4)
plot(InputRspRange,InputRspDist)
xlim([min(InputRspRange)-1,max(InputRspRange)+1])
subplot(2,2,1)
plot(OutputRspDist,OutputRspRange)
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




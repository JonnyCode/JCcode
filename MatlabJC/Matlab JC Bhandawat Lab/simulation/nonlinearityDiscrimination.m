% simulation of nonlinearity;s effect on discriminability
% JC 9/27/12
clear all

Xdelta = .1 ;
Xmax = 50 ;

g1_Xmean = 25 ;
g1_Xstd = 4 ;

g2_Xmean = 27 ;
g2_Xstd = 2 ;

Xrange = [0:Xdelta:50] ;
g1 = exp(-((Xrange-g1_Xmean).^2)/(2*g1_Xstd^2)) ; % gaussian
g2 = exp(-((Xrange-g2_Xmean).^2)/(2*g2_Xstd^2)) ; % gaussian
g1= g1/sum(g1) ;
g2 = g2/sum(g2) ;
gin = [g1;g2] ;

% NL = g2./(g1+g2) ;
% NL = NL*Xmax ;
NL = cumsum(g1) ;

% NL
[OutputDist,OutputDistX] = pdfTransform(g1,Xrange,NL) ;


NLunique = unique(NL) ;
NLdiffSmooth = smooth([diff(NLunique(1:2)),diff([NLunique])],3) ;
for a=1:length(NLunique) ;
    g1out(a) = sum(g1.*(NL==NLunique(a)))/NLdiffSmooth(a) ;
    g2out(a) = sum(g2.*(NL==NLunique(a)))/NLdiffSmooth(a) ;
end
gout = [g1out; g2out] ;    

figure
plot(NLunique,g1out) ;













% suming weighted distributions
g1inSum = sum(Xrange.*g1) 
g2inSum = sum(Xrange.*g2) 

g1outSum = sum(NLunique.*g1out) 
g2outSum = sum(NLunique.*g2out) 

% discriminability
Dgin = sum(max(g,[],1))/sum(sum(g)) ;
Dgout = sum(max(gout,[],1))/sum(sum(gout)) ;

% SNR
g1inMean = sum(g1.*Xrange) ;
g1inVar = sum(g1.*(Xrange-g1inMean).^2)

g1Var = var(g1out) ;
g2Mean = mean(g2out)

close all

figure
plot(Xrange,g1)
hold on
plotyy(Xrange,g2,Xrange,NL)
plot(Xrange,g2,'r')

figure
plot(NLunique,g1out)
hold on
plot(NLunique,g2out,'r')






%output distributions 
OutputRspDist = nan(2,length(Xrange)) ;
OutputRspRange = NL ;
a=1 ;
b=1 ;
while a<=length(Xrange) ;
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
        OutputRspDist(:,b) = (g(:,a:a-1+i)*weightingFnct'*Xdelta/NL(a)) ; 
        a = a + i  ;
        
    elseif a~=length(Xrange) ; 
        i = find(NL(a+1:end)~=NL(a),1,'first') ;
        if isempty(i) ;
            i = length(NL)-a+1 ;
        end
        weightingFnct = [1,2*ones(1,i-1),1] ; % NEEDS TO BE FIXED
        OutputRspDist(:,b) = (g(:,a-1:a-1+i)*weightingFnct'*Xdelta/(NL(a)-NL(a-1)))-OutputRspDist(:,b-1) ;
        a = a + i  ;
    end
    
    b=b+1 ;
end

OutputRspDist_sum = sum(OutputRspDist,1) ;


% output distributions
NL = NL + [10^-10:10^-10:length(NL)*10^-10] ; % this is a hack and should be removed once issue is resolved properly

OutputRspDist = nan(2,length(InputRspRange)) ;
for a=1:length(InputRspRange) ;
    for b=1:2 ;
        if a==1 ;        
            OutputRspDist(b,a) = InputRspDist(b,a)*Xdelta/NL(a) ;
            OutputRspRange(a) = NL(a) ;
            
        else ;
            OutputRspDist(b,a) = ((InputRspDist(b,a-1)+InputRspDist(b,a))*Xdelta/abs(NL(a-1)-NL(a)))-OutputRspDist(b,a-1) ;
            OutputRspRange(a) = NL(a) ;
        end
    end
end




% simulation of anticorrelation bias as a function of sample number

figure
for trial = 1:100 ;

SampleNumRange = [2:1:50] ;


num=0 ;
for n=SampleNumRange ;
num = num+1 ;
    
x=normrnd(100,60,1,n) ;

index = 1:length(x) ;
for a=1:length(x) ;
    R(a) = x(a)-mean(x) ;
    R2(a) = x(a)-mean(x(index(index~=a))) ;
end

rnd=0 ;
for a=1:length(x)-1 ;
    for b = a+1:length(x) ;
        rnd=rnd+1 ;
        c(rnd) = (R(a)*R(b))/(R(a)^2+R(b)^2) ;
        c2(rnd) = (R2(a)*R2(b))/(R2(a)^2+R2(b)^2) ;
    end
end

cc(trial,num) = mean(c) ;
cc2(trial,num) = mean(c2) ;

end
clearvars -except cc cc2 SampleNumRange trial

plot(SampleNumRange,cc(trial,:))
hold on
plot(SampleNumRange,cc2(trial,:),'k')

end

plot(SampleNumRange,mean(cc),'r')
plot(SampleNumRange,mean(cc2),'g--')

%

SampleNumRange = [4:1:50] ;

for trial = 1:length(SampleNumRange) ;

    for a=1:SampleNumRange(trial) ;
        exc(a,:) = normrnd(100,60,1,100000) ;
        inh(a,:) = exc(a,:) ;
    end

    shuffexc = circshift(exc,[2,0]) ;
    shuffinh = circshift(inh,[3,0]) ;
    
    index = 1:SampleNumRange(trial) ;
    for a=1:SampleNumRange(trial) ;
        Re(a,:) = exc(a,:) -   ;
        Ri(a,:) = inh(a,:) - shuffinh(a,:)  ;
    end
    
    shufRi = circshift(Ri,[1,0]) ;
    
    for a=1:SampleNumRange(trial) ;
        ccShuff(a,:) = xcov(Re(a,:),shufRi(a,:),'coef') ;
        ccSim(a,:) = xcov(Re(a,:),Ri(a,:),'coef') ;
    end
    
    ccMeanShuff(trial,:) = mean(ccShuff) ;
    corrcoefShuff(trial) = ccMeanShuff(trial,length(Re)) ;
    
    ccMeanSim(trial,:) = mean(ccSim) ;
    corrcoefSim(trial) = ccMeanSim(trial,length(Re)) ;
end

SampleNumRange = [2:1:50] ;

for trial = 1:length(SampleNumRange) ;

    for a=1:SampleNumRange(trial) ;
        exc(a,:) = normrnd(100,60,1,100000) ;
        inh(a,:) = exc(a,:) ;
    end

    exc_mean = mean(exc) ;
    inh_mean = mean(inh) ;
    
    for a=1:SampleNumRange(trial) ;
        cc(a,:) = xcov(exc(a,:),inh(a,:)) ;
        
        ac_exc(a,:) = xcov(exc(a,:)) ;
        ac_inh(a,:) = xcov(inh(a,:)) ;
    end
    
    cc_mean = mean(cc) ;
    
    ac_exc_mean = mean(ac_exc) ;
    ac_inh_mean = mean(ac_inh) ;
    vartotal = max(sqrt(ac_exc_mean).*sqrt(ac_inh_mean)) ;
    
    Meancc(trial,:) = xcov(exc_mean,inh_mean) ;
    
    cc_shuffcorrected(trial,:) = (cc_mean - Meancc(trial,:))./(vartotal - Meancc(trial,:)) ;
    
end    
    





pnts = 100000
anticorr = 1 ; % 1 if signals are correlated, negative 1 if signals are anticorrelated


signal = normrnd(50,100,1,pnts) ;
signal = lowPassFilter(signal, 10000, 50) ; % low pass filter

shift = 0:10:100 ;

for a =1:length(shift);

    for trial = 1:10 ;
        
        x=normrnd(50,5,1,pnts) ; % mean,std,m,n
        y=normrnd(50,5,1,pnts) ;
        s=normrnd(50,100,1,pnts) ;
        s = lowPassFilter(s, 10000, 50) ; % low pass filter
    
        xx(trial,:)=x+s ;
        yy(trial,:)=(y+circshift(s,[0,shift(a)]))*anticorr ;
        
        zz(trial,:)=xx(trial,:)-yy(trial,:)+signal ;

    end
    
    for trial = 1:10 ;
        xx_Res(trial,:) = mean(xx)-xx(trial,:) ;
        yy_Res(trial,:) = mean(yy)-yy(trial,:) ;
        cc(trial,:) = xcorr(xx_Res(trial,:),yy_Res(trial,:),'coef') ;
    end
    ccMean(a) = mean(cc(:,pnts)) ;  
    
    snrMean(a) = mean(mean(zz)./std(zz)) ;
   
    k = mean(mean(zz)./sqrt(var(xx)+var(yy))) ;
    snrPred(a) = k/sqrt(1-ccMean(a)) ;

end

figure
plot(ccMean,snrMean,'b-') ;
hold on
plot(ccMean,snrPred,'r--') ;

pnts = 100000
anticorr = -1 ; % 1 if signals are correlated, negative 1 if signals are anticorrelated

x=normrnd(0,5,1,pnts) ; % mean,std,m,n
y=normrnd(0,5,1,pnts) ;
s=normrnd(0,100,1,pnts) ;

s = lowPassFilter(s, 10000, 50) ; % low pass filter

round = 0 ;

for shift = 0:10:1000 ;

round = round+1 ;
    
xx(round,:)=x+s ;
yy(round,:)=(y+circshift(s,[0,shift]))*anticorr ;
zz(round,:)=xx(round,:)+yy(round,:) ;

c(round,:) = xcorr(xx(round,:),yy(round,:)) ;

varxx(round) = var(xx(round,:)) ;
varyy(round) = var(yy(round,:)) ;
vars = var(s) ;

varzzPred(round) = varxx(round) + varyy(round) + (2*c(round,pnts))/pnts ;

varzz(round) = var(zz(round,:)) ;

end

figure
plot(0:10:1000,varzz)
hold on
plot(0:10:1000,varzzPred,'r--')
title('simplesim')
legend('var z', 'analytical')
xlabel('time shift')
ylabel('variance')

figure
plot(2*c(:,pnts)/pnts,varzzPred,'b*') ;






pnts = 100000
anticorr = -1 ; % 1 if signals are correlated, negative 1 if signals are anticorrelated
stds(1) = 5 ;

for level=1:10 ;

x=normrnd(0,5,1,pnts) ; % mean,std,m,n
y=normrnd(0,5,1,pnts) ;
s=normrnd(0,stds(levels),1,pnts) ;

s = lowPassFilter(s, 10000, 50) ; % low pass filter

for shift = 10 ;
   
xx(level,:)=x+s ;
yy(level,:)=(y+circshift(s,[0,shift]))*anticorr ;
zz(level,:)=xx(level,:)+yy(level,:) ;

c(level,:) = xcorr(xx(level,:),yy(level,:)) ;

varxx(level) = var(xx(level,:)) ;
varyy(level) = var(yy(level,:)) ;
vars = var(s) ;

varzzPred(round) = varxx(level) + varyy(level) + (2*c(level,pnts))/pnts ;

varzz(level) = var(zz(level,:)) ;


end
stds(level+1)= sqrt(varzz(level)) ;

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






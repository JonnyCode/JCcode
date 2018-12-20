% sigmoid against sigmoid

Xmax = 20 ;
Xdelta = .1 ;

g1_Xmean = 6 ;
g1_Xstd = 2 ;

g2_Xmean = 6 ;
g2_Xstd = 3 ;

g3_Xmean = 5 ;
g3_Xstd = 1 ;

Xrange = [0:Xdelta:Xmax] ;

g1 = exp(-((Xrange-g1_Xmean).^2)/(2*g1_Xstd^2)) ; % gaussian
g2 = exp(-((Xrange-g2_Xmean).^2)/(2*g2_Xstd^2)) ; % gaussian
g3 = exp(-((Xrange-g3_Xmean).^2)/(2*g3_Xstd^2)) ; % gaussian

s1 = cumsum(g1)/sum(g1) ;
s2 = cumsum(g2)/sum(g2) ;
s3 = cumsum(g3)/sum(g3) ;

figure
plot(Xrange,s1)
hold on
plot(Xrange,s2,'r')
plot(Xrange,s3,'g')

figure
plot(s1,s2,'k')
hold on
plot(s1,s3,'k--')


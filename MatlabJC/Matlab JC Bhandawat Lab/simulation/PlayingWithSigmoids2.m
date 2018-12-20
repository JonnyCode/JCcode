% sigmoid against sigmoid

Xmax = 20 ;
Xdelta = .1 ;

g1_Xmean = [7,7,7.5,8.5,10] ;
g1_Xstd = [2,2,2,2,2] ;
g1_amp = [3,3,3,2.5,2.4] ;

g2_Xmean = [2,4,6,8,10] ;
g2_Xstd = [3.75,3,2.5,2.5,2] ;
g2_amp = [1,1,1,.8,.6] ;

Xrange = [0:Xdelta:Xmax] ;

for a=1:length(g1_Xmean) ;

    g1 = exp(-((Xrange-g1_Xmean(a)).^2)/(2*g1_Xstd(a)^2)) ; % gaussian
    g2 = exp(-((Xrange-g2_Xmean(a)).^2)/(2*g2_Xstd(a)^2)) ; % gaussian

    s1(a,:) = g1_amp(a)*cumsum(g1)/sum(g1) ;
    s2(a,:) = g2_amp(a)*cumsum(g2)/sum(g2) ;
    
    plot(s1(a,:),s2(a,:))
    hold on
end

figure
plot(Xrange,s1)
hold on
plot(Xrange,s2,'--')

figure
for a=1:length(g1_Xmean) ;
    plot(s1(a,:),s2(a,:),'k','LineWidth',a)
    hold on
end


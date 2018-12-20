% Gaussian Cartoon for figure 5 of pair ds paper

% JC 10/26/11
X = [0:.01:10] ;

std1=1 ;
std2=1.2 ;

mean1 = 3 ;
mean2 = 6 ;

cov1=(std1-.3)^2 ;

h=error_ellipse([std1^2,0;0,std1^2],[mean1,mean1],'conf',.95) ; % confidence elipse
ellipse1X = get(h,'XData') ;
ellipse1Y = get(h,'YData') ;

h=error_ellipse([std1^2,0;0,std1^2],[mean2,mean2],'conf',.95) ; % confidence elipse
ellipse2X = get(h,'XData') ;
ellipse2Y = get(h,'YData') ;

h=error_ellipse([std2^2,0;0,std2^2],[mean1,mean1],'conf',.95) ; % confidence elipse
ellipse3X = get(h,'XData') ;
ellipse3Y = get(h,'YData') ;

h=error_ellipse([std2^2,0;0,std2^2],[mean2,mean2],'conf',.95) ; % confidence elipse
ellipse4X = get(h,'XData') ;
ellipse4Y = get(h,'YData') ;


h=error_ellipse([std1^2,cov1;cov1,std1^2],[mean1,mean1],'conf',.95) ; % confidence elipse
ellipse5X = get(h,'XData') ;
ellipse5Y = get(h,'YData') ;

h=error_ellipse([std1^2,cov1;cov1,std1^2],[mean2,mean2],'conf',.95) ; % confidence elipse
ellipse6X = get(h,'XData') ;
ellipse6Y = get(h,'YData') ;


h=error_ellipse([std1^2,0;0,std1^2],[mean1,mean2],'conf',.95) ; % confidence elipse
ellipse7X = get(h,'XData') ;
ellipse7Y = get(h,'YData') ;

h=error_ellipse([std1^2,0;0,std1^2],[mean2,mean1],'conf',.95) ; % confidence elipse
ellipse8X = get(h,'XData') ;
ellipse8Y = get(h,'YData') ;

h=error_ellipse([std1^2,cov1;cov1,std1^2],[mean1,mean2],'conf',.95) ; % confidence elipse
ellipse9X = get(h,'XData') ;
ellipse9Y = get(h,'YData') ;

h=error_ellipse([std1^2,cov1;cov1,std1^2],[mean2,mean1],'conf',.95) ; % confidence elipse
ellipse10X = get(h,'XData') ;
ellipse10Y = get(h,'YData') ;


%figure
figure
plot(ellipse1X,ellipse1Y)
hold on
plot(ellipse2X,ellipse2Y)
plot(ellipse3X,ellipse3Y,'--')
plot(ellipse4X,ellipse4Y,'--')

figure
plot(ellipse1X,ellipse1Y,'--')
hold on
plot(ellipse2X,ellipse2Y,'--')
plot(ellipse5X,ellipse5Y)
plot(ellipse6X,ellipse6Y)

figure
plot(ellipse7X,ellipse7Y,'--')
hold on
plot(ellipse8X,ellipse8Y,'--')
plot(ellipse9X,ellipse9Y)
plot(ellipse10X,ellipse10Y)

% for igor

identifier = ['ellipse1X'] ;
ForIgor.(identifier) = ellipse1X ; 

identifier = ['ellipse1Y'] ;
ForIgor.(identifier) = ellipse1Y ; 

identifier = ['ellipse2X'] ;
ForIgor.(identifier) = ellipse2X ; 

identifier = ['ellipse2Y'] ;
ForIgor.(identifier) = ellipse2Y ; 

identifier = ['ellipse3X'] ;
ForIgor.(identifier) = ellipse3X ; 

identifier = ['ellipse3Y'] ;
ForIgor.(identifier) = ellipse3Y ; 

identifier = ['ellipse4X'] ;
ForIgor.(identifier) = ellipse4X ; 

identifier = ['ellipse4Y'] ;
ForIgor.(identifier) = ellipse4Y ; 

identifier = ['ellipse5X'] ;
ForIgor.(identifier) = ellipse5X ; 

identifier = ['ellipse5Y'] ;
ForIgor.(identifier) = ellipse5Y ; 

identifier = ['ellipse6X'] ;
ForIgor.(identifier) = ellipse6X ; 

identifier = ['ellipse6Y'] ;
ForIgor.(identifier) = ellipse6Y ; 

identifier = ['ellipse7X'] ;
ForIgor.(identifier) = ellipse7X ; 

identifier = ['ellipse7Y'] ;
ForIgor.(identifier) = ellipse7Y ; 

identifier = ['ellipse8X'] ;
ForIgor.(identifier) = ellipse8X ; 

identifier = ['ellipse8Y'] ;
ForIgor.(identifier) = ellipse8Y ; 

identifier = ['ellipse9X'] ;
ForIgor.(identifier) = ellipse9X ; 

identifier = ['ellipse9Y'] ;
ForIgor.(identifier) = ellipse9Y ; 

identifier = ['ellipse10X'] ;
ForIgor.(identifier) = ellipse10X ; 

identifier = ['ellipse10Y'] ;
ForIgor.(identifier) = ellipse10Y ; 






% anti correlated traces that are correlated by a discreet event passing
% threshold
time = [.0001:.0001:1] ; 

trace = normrnd(0,1,1,length(time)) ;
trace = lowPassFilter(trace,10000,100) ;
[m,im] = max(trace(2000:8000)) ;

correlatingBump = zeros(1,length(time)) ;
correlatingBump(im+2000) = 2 ;
correlatingBump = lowPassFilter(correlatingBump,10000,100) ;

correlatingBump1 = correlatingBump + 1 ;
correlatingBump2 = ((correlatingBump/max(correlatingBump))*-(max(correlatingBump)+2))+1 ;

trace1 = trace.*correlatingBump1 ;
trace2 = -trace.*correlatingBump2 ;

% correlated traces that are uncorrelated by a two antiphase sinwaves

trace3 = trace+sind(50000*time).*.05 ;
trace4 = trace+sind(50000*time+181).*.05 ;

% peaksi1 = localMaxFinder(trace1) ;
% peaksi2 = localMaxFinder(trace2) ;
% 
% peaks = sort([trace1(peaksi1),trace2(peaksi2)]) ;
% Threshold = peaks(end-2) ;

figure
plot(time,trace1)
hold on
plot(time,trace2,'r')

figure
plot(time,trace3)
hold on
plot(time,trace4,'r')


% for Igor

identifier = ['trace1'] ;
ForIgor.(identifier) = trace1 ; 

identifier = ['trace2'] ;
ForIgor.(identifier) = trace2 ; 

identifier = ['trace3'] ;
ForIgor.(identifier) = trace3 ; 

identifier = ['trace4'] ;
ForIgor.(identifier) = trace4 ; 








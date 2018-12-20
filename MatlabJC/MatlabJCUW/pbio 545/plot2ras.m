function plot2ras(row,intervals,nsteps)

%  Support function for GamSpikeTrainTutorial.m
% 
%  Just so we don't have to copy this code a million times.
%  6/14/00 GDLH
%  1/28/03 modified slightly, mns

x = [-5:.1:60];
subplot(3,2,row*2-1);
[n,bins] = hist(intervals,x);
n = n/(sum(n)*.1);
bar(x(1:end-1),n(1:end-1));
set(gca,'XLim',[min(x) max(x)]);
set(gca,'Ylim',[0 .18])
ylabel('Probability')
xlabel('Interval')
if (nsteps == 1)
	text(20,.15,[num2str(nsteps),' step to threshold']);
else
	text(20,.15,[num2str(nsteps),' steps to threshold']);
end

subplot(3,2,row*2);
arrivals = cumsum(intervals);
plot1ras(arrivals(arrivals<500))
set(gca,'Xlim',[0 500]);


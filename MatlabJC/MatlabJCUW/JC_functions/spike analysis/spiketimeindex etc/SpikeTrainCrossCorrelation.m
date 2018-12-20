function [AvgCrossCorrelation,xvalues,sigma3]=SpikeTrainCrossCorrelation(CellParameters)


% for cell attached data

clear AvgCrossCorrelation CrossCorrelationMatrix SpikeTrain comparisoncount i j k l m n

for i=1:length(CellParameters);                   % for each of the trials at this mean light intensity....
    SpikeTrain(i,:)=zeros(1,75000);               % ...make a new matrix composed of zeros
    for j=1:length(CellParameters{i});            % for each of the spikes in this particular trial...
        a=CellParameters{i}(j);                   % ....at what time did the spike occur
        SpikeTrain(i,ceil(a))=1;                  % at the time the spike occurred (ROUNDED TO THE NEAREST MS), replace the 0 with a 1
        clear a;                                                        
    end
end


[m,n]=size(SpikeTrain);
comparisoncount=1;
for k=1:m-1;
    for l=k+1:m;
        CrossCorrelationMatrix(comparisoncount,:)=xcorr(SpikeTrain(k,:), SpikeTrain(l,:), 50, 'coeff');            % make a new matrix with cross-correlations from all possible comparisons of the spike trains
        comparisoncount=comparisoncount+1;
    end
end

AvgCrossCorrelation=mean(CrossCorrelationMatrix,1);
[c,d]=size(AvgCrossCorrelation);
 
xvalues=-1*((d-1)/2):1:((d-1)/2);

coef=[1 1 15];
fitcoef=nlinfit(xvalues,AvgCrossCorrelation,'gaussian',coef);
fit=gaussian(fitcoef,xvalues);

plot(xvalues,AvgCrossCorrelation,'bo')
hold on
plot(xvalues,fit,'k')

sigma3=abs(round(fitcoef(3)*3));



% % for whole cell data
%       
% [GoodEpochData, Offset]=GetGoodEpochData(CellInfo, EpochCondition(3));
% [m,n]=size(GoodEpochData);
% comparisoncount=1;
% for k=1:m-1;
%     for l=k+1:m;
%         CrossCorrelationMatrix(comparisoncount,:)=xcorr(GoodEpochData(k,:), GoodEpochData(l,:), 500, 'coeff');            % make a new matrix with cross-correlations from all possible comparisons of the spike trains
%         comparisoncount=comparisoncount+1;
%     end
% end
% 
% AvgCrossCorrelation=mean(CrossCorrelationMatrix,1)

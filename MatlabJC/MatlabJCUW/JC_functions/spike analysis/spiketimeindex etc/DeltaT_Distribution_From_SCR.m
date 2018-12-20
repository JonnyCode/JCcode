function [DeltaT, tli_spike, tlj_spike, Percent_Pairs_Quantified, DeltaT_Histogram,X_Values,DeltaT_CumProb] = DeltaT_Distribution_From_SCR(SpikeTimeIndex,cost)



%
% determines the distribution of delta t values for each
% valid pairwise spike train comparison (from a cell array that contains 
% x number of spike train spike times)
%


numbertrials=length(SpikeTimeIndex);
comparisoncount=1;                                      % # of spike train comparisons performed

for j=1:numbertrials-1;                                      % for all of the spike trains except the last...
    for k=j+1:numbertrials;                                    % for the spike train # greater than j
        a=length(SpikeTimeIndex{j});                      % determine the # of spikes in spike train 1 
        b=length(SpikeTimeIndex{k});                        % determine the # of spikes in spike train 2 
        [d,scr]=spkd_with_scr(SpikeTimeIndex{j},SpikeTimeIndex{k},cost);                    % compute scr matrix
        [DeltaT{comparisoncount}, tli_spike{comparisoncount}, tlj_spike{comparisoncount}]=DeltaT_From_SCR(scr,cost,SpikeTimeIndex{j},SpikeTimeIndex{k});
        TotalNumberOfSpikes(comparisoncount)=a+b;                            % compute the total number of spikes in the two spike trains being compared
        Percent_Pairs_Quantified(comparisoncount)=length(DeltaT{comparisoncount}).*2/TotalNumberOfSpikes(comparisoncount);
        comparisoncount=comparisoncount+1;
        clear b;
    end
    clear a;
end

DeltaT_Histogram=zeros(comparisoncount-1,1000);

for a=1:comparisoncount-1                                                                                                                        % for each comparison
        for X_Values=0:1:999                                                                                                                            % for each DeltaT interval less than 100ms 
        temp1=find(DeltaT{a}>(((X_Values)/10)-.001) & DeltaT{a}<(((X_Values)/10)+.001));                % find the instances where spikes are separated by that interval
        temp2=length(temp1);                                                                                                                                 % count the number of instances
        DeltaT_Histogram(a,X_Values+1)=temp2;                                                                                               % place that number in the appropriate slot
        clear temp*
    end
end

X_Values=0:0.1:99.9;
Summed_DeltaT_Histogram=sum(DeltaT_Histogram,1);
DeltaT_CumProb=cumsum(Summed_DeltaT_Histogram)./max(cumsum(Summed_DeltaT_Histogram));


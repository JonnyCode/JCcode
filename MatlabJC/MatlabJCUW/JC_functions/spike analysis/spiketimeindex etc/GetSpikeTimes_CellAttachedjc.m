function [SpikeTimeIndex]=GetSpikeTimes_CellAttachedjc(CellInfo, Condition, pre, post, threshold)
%
% Obtain time for all events with an ampltidue >= a set threshold,
% excluding those that occur before and after stimulus; when events
% occur within 1000/samplinginterval points(e.g., 10 points, or 1ms,
% if sampling interval=100 us), only the first point above threshold 
% is retained. * this last point is changed in revision (see below) *
%
% GJM 11/04
%    refined with help of GF and FR 12/04
%       
% GJM 1/05
%   * added second criteria for spike determination - there must be a positive deflection >threshold 1-3 points
%   after negative deflection (that is itself more negative than negative threshold) 
%
%   * made into function that can be called form command line
%
%   * corrected two major errors:
%               1) routine throwing out spikes in which more than one point in the spike waveform was above threshold
%               2) routine not throwing out putative spikes that were separated by < refractory period    
%       *as a result, all spike distance measures computed before 1/25/05 must be calculated again*
%
% GJM 12/05
%   * spikes now timed at their peak, rather than at the first point they cross the negative threshold 
%           (should help avoid differences in timing due to differences in spike height)
%                   all data prior to 120105 needs to be reanalyzed
%

% clear TempData TempSpikeTime SpikeTimeIndex  pre firststimpt post laststimpt tempnumberspikes f samplinginterval epochs totalepochpoints

LowPassFilteredData = ApplyHiFrequencyCutoff(CellInfo, Condition, 20);   % remove low frequency noise
[GoodEpochData, Offset]=GetGoodEpochData(CellInfo, LowPassFilteredData);          % remove trials flagged during parsing

[epochs,totalepochpoints]=size(GoodEpochData);

% pre=0
%pre=FindSearchPara(Condition, 'PrePoints');         % # of values collected before stimulus
firststimpt=pre+1;                                          % first value during stimulus

% post=0
%post=FindSearchPara(Condition, 'TailPoints');       % # of values collected after stimulus
laststimpt=totalepochpoints-post;                                   % last value during stimulus
    
samplinginterval=FindSearchPara(Condition,'SampInterv');        % sampling interval (in microseconds)
timefactor=(1/samplinginterval)*1000;                                    % scales points into ms timeframe (if sampling interval is 100us, there are 10 points per ms)
samplingintervalinms=samplinginterval/1000   ;                                                                     

% now identify negative deflections < -threshold that are followed by quickly by a
% positive deflection > threshold

                                                                 
for o = 1:epochs                                                     % for # of identical stimulus repeats....
    TempData(o,1:laststimpt-pre) = GoodEpochData(o,firststimpt:laststimpt);     % move only those points collected during stimulus into new matrix
    TempSpikeTime = find(TempData(o,:)<-1*threshold)';                          % find points collected during stimulus (in a single epoch) where value is < -1*threshold
    [tempnumberspikes,f]=size(TempSpikeTime);                                   % find number of suprathreshold points during this epoch (ie, i=1, i=2, etc.)
    count=1;                                                                        
    for t=1:tempnumberspikes                                          % for each point (except the last) that passes threshold in this epoch ....
        if (TempSpikeTime(t)<(laststimpt-pre)-5) & (TempData(o,1+TempSpikeTime(t))>threshold) ;                    % ...if the next point is above the positive threshold....
           TempSpikeTimeIndex1{o}(count)=TempSpikeTime(t);        % ...keep the point - index it according to how many points have already been retained
           count=count+1;                                               % update counter 
        elseif (TempSpikeTime(t)<(laststimpt-pre)-5) & (TempData(o,2+TempSpikeTime(t))>threshold) ;                     %...if the point after the next point is above the positive threshold....
           TempSpikeTimeIndex1{o}(count)=TempSpikeTime(t);        % ...keep the point - index it according to how many points have already been retained
           count=count+1; 
        elseif (TempSpikeTime(t)<(laststimpt-pre)-5) & (TempData(o,3+TempSpikeTime(t))>threshold) ; 
           TempSpikeTimeIndex1{o}(count)=TempSpikeTime(t);        
           count=count+1; 
        elseif (TempSpikeTime(t)<(laststimpt-pre)-5) & (TempData(o,4+TempSpikeTime(t))>threshold) ; 
           TempSpikeTimeIndex1{o}(count)=TempSpikeTime(t);        
           count=count+1;    
        elseif (TempSpikeTime(t)<(laststimpt-pre)-5) & (TempData(o,5+TempSpikeTime(t))>threshold) ; 
           TempSpikeTimeIndex1{o}(count)=TempSpikeTime(t);        
           count=count+1;    
        end
    end
     clear TempSpikeTime
end
  

%
% now save only a single spike timer per spike event; don't want to save two
% time points that simply represent the same spike (at two points which are
% both below threshold).
% 


for b=1:length(TempSpikeTimeIndex1)                      % for each epoch 
    count=1;
    Diff_Times=diff(TempSpikeTimeIndex1{b});            % calculate the difference in time between points below the negative threshold
    for stepper=1:length(Diff_Times);                    % for each of those points
         if Diff_Times(stepper)>1.5;                      % if the difference between one point and the next is more than .15ms
             transition(count)=stepper;                     % consider it the last negative point of a putative - make a list of these transition times
             count=count+1;
         end
    end
    for spike=1                                                                                         % for the first putative spike
        [temp1(spike),temp2(spike)]=min(TempData(b,TempSpikeTimeIndex1{b}(1:transition(spike))));     % determine the time at which the spike reaches its negative peak
        TempSpikeTimeIndex2{b}(spike)=TempSpikeTimeIndex1{b}(temp2(spike))/10;                               % put the time, converted back to ms, in the index  
    end
    for spike=2:length(transition)                                                                       % for each putative spike except the last
        [temp1(spike),temp2(spike)]=min(TempData(b,TempSpikeTimeIndex1{b}(transition(spike-1)+1:transition(spike))));
        TempSpikeTimeIndex2{b}(spike)=TempSpikeTimeIndex1{b}(transition(spike-1)+temp2(spike))/10;
    end
    for spike=length(transition)+1
        [temp1(spike),temp2(spike)]=min(TempData(b,TempSpikeTimeIndex1{b}(transition(spike-1)+1:end)));
        TempSpikeTimeIndex2{b}(spike)=TempSpikeTimeIndex1{b}(transition(spike-1)+temp2(spike))/10;
    end
    clear transition stepper Diff_Times
end
    
      
% now check that spikes are seperated by more than the refractory period
% (>1ms)


[f,g]=size(TempSpikeTimeIndex2);
for y=1:g;
    [h]=length(TempSpikeTimeIndex2{y});
    count=1;
    for z=count:h
        if z==1; 
            SpikeTimeIndex{y}(count)=TempSpikeTimeIndex2{y}(z);
            count=count+1;
        else z~=1 & z<=h;
            if TempSpikeTimeIndex2{y}(z)-TempSpikeTimeIndex2{y}(z-1)>=10*samplingintervalinms;
            SpikeTimeIndex{y}(count)=TempSpikeTimeIndex2{y}(z);
            count=count+1;
            end
        end
    end
end


    
    















%--------------------------------




% EARLY VERSION
%   Contains a very bad error - relevant spikes are thrown out (if two points in the same spike are above/below threshold) 
%
%
% clear TempData TempSpikeTime SpikeTimeIndex  pre firststimpt post laststimpt tempnumberspikes f samplinginterval epochs totalepochpoints
% 
% LowPassFilteredData = ApplyHiFrequencyCutoff(CellInfo, EpochCondition(1), 20);   % remove low frequency noise
% [GoodEpochData, Offset]=GetGoodEpochData(CellInfo, LowPassFilteredData);          % remove trials flagged during parsing
% 
% [epochs,totalepochpoints]=size(GoodEpochData);
% 
% pre=FindSearchPara(EpochCondition(1), 'PrePoints');         % # of values collected before stimulus
% firststimpt=pre+1;                                          % first value during stimulus
% 
% post=FindSearchPara(EpochCondition(1), 'TailPoints');       % # of values collected after stimulus
% laststimpt=totalepochpoints-post;                                   % last value during stimulus
%     
% samplinginterval=FindSearchPara(EpochCondition(1),'SampInterv');        % sampling interval (in microseconds)
% timefactor=(1/samplinginterval)*1000                                    % scales points into ms timeframe (if sampling interval is 100us, there are 10 points per ms)
% 
% for o = 1:epochs                                                      % for # of identical stimulus repeats....
%     TempData(o,1:laststimpt-pre) = GoodEpochData(o,firststimpt:laststimpt);     % move only those points collected during stimulus into new matrix
%     TempSpikeTime = find(TempData(o,:)>10)';                                    % find points collected during stimulus (in a single epoch) where value is >= threshold
%     [tempnumberspikes,f]=size(TempSpikeTime);                                   % find number of suprathreshold points during this epoch (ie, i=1, i=2, etc.)
%     count=1;                                                                        
%     for t=1:tempnumberspikes-1                                          % for each point above threshold in this epoch (except the last)....
%         if TempSpikeTime(t+1)-TempSpikeTime(t)>1000/samplinginterval;   % ...if the next point above threshold is more than 1 ms away....
%            SpikeTimeIndex{o}(count)=TempSpikeTime(t)/timefactor;        % ...keep the point (and convert its time to ms) - index it according to how many points have already been retained
%            count=count+1;                                               % update counter if suprathreshold points 1 & 2 were more than 1ms apart
%         end
%     end
%     clear TempSpikeTime TempData
% end
% 
% % then need to pull data from cell, perform computation of choice, and put result into matrix
% % for example, compute average 'victor spike distance' for a given epoch condition by making 
% % all valid pairwise comparisons. 
% 
% clear DMatrix AvgCostMatrix cost costindex
% comparisoncount=1                                       % # of spike train comparisons performed
% for j=1:epochs-1;                                       % for all of the spike trains except the last...
%     for k=j+1:epochs;                                   % 
%         for costindex=1:50;                          
%             cost(costindex)=costindex/2;
%             DMatrix(comparisoncount,costindex)=spkd(SpikeTimeIndex{j}, SpikeTimeIndex{k}, (cost(costindex)));   %compute spike distance for a given cost; place value into appropriate slot in DMatrix
%         end
%         comparisoncount=comparisoncount+1
%     end
% end
% 
% AvgCostMatrix=mean(DMatrix,1)       % compute average spike distance across all pairwise comparisons for a given cost 
% 
% 
% CellParameters.CostMatrix=DMatrix
% CellParameters.Cost=cost
% CellParameters.AvgCostMatrix=AvgCostMatrix
% 
% % add pop-up prompt here to save



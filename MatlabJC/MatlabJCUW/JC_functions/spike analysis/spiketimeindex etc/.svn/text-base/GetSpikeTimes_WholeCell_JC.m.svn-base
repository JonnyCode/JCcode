function [SpikeTimeIndex,totalepochpoints]=GetSpikeTimes_WholeCell_JC(CellInfo,Condition,threshold)
%
% Obtain time for all events with an ampltidue >= a set threshold,
% excluding those that occur before and after stimulus; when events
% occur within 1000/samplinginterval points(e.g., 10 points, or 1ms,
% if sampling interval=100 us), only the first point above threshold 
% is retained.
%
%
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
%   * spikes now timed at their peak - should help avoid differences in timing due to differences in spike height)
%                   all data prior to 120105 needs to be reanalyzed


clear TempData TempSpikeTime SpikeTimeIndex  pre firststimpt post laststimpt tempnumberspikes f samplinginterval epochs totalepochpoints

[GoodEpochData_prelim, Offset]=GetGoodEpochData(CellInfo, Condition);          % remove trials flagged during parsing

[epochs,totalepochpoints]=size(GoodEpochData_prelim);

for trials=1:epochs
     GoodEpochData(trials,:)=GoodEpochData_prelim(trials,:)+Offset(trials);
end

pre=0
%pre=FindSearchPara(Condition, 'PrePoints');         % # of values collected before stimulus
firststimpt=pre+1;                                          % first value during stimulus

post=0
%post=FindSearchPara(Condition, 'TailPoints');       % # of values collected after stimulus
laststimpt=totalepochpoints-post;                                   % last value during stimulus
    
samplinginterval=FindSearchPara(Condition,'SampInterv');        % sampling interval (in microseconds)
timefactor=(1/samplinginterval)*1000;                                    % scales points into ms timeframe (if sampling interval is 100us, there are 10 points per ms)
samplingintervalinms=samplinginterval/1000   ;                                                                     

% now identify negative deflections < -threshold that are followed by quickly by a
% positive deflection > threshold

                                                                 
for o = 1:epochs                                                     % for # of identical stimulus repeats....
    TempData(o,1:laststimpt-pre) = GoodEpochData(o,firststimpt:laststimpt);     % move only those points collected during stimulus into new matrix
    TempSpikeTime = find(TempData(o,:)>threshold)';                          % find points collected during stimulus (in a single epoch) where value is > threshold
    [tempnumberspikes,f]=size(TempSpikeTime);                                   % find number of suprathreshold points during this epoch (ie, i=1, i=2, etc.)
    count=1;                                                                        
    for t=1:tempnumberspikes                                          % for each point (except the last) that passes threshold in this epoch ....
        if (TempSpikeTime(t)<(laststimpt-pre)-5) & (TempData(o,1+TempSpikeTime(t))<threshold) ;                    % ...if the next point is below the threshold....
           TempSpikeTimeIndex1{o}(count)=TempSpikeTime(t);        % ...keep the point - index it according to how many points have already been retained
           count=count+1;                                               % update counter 
        elseif (TempSpikeTime(t)<(laststimpt-pre)-5) & (TempData(o,2+TempSpikeTime(t))<threshold) ;                     %...if the point after the next point is below the threshold....
           TempSpikeTimeIndex1{o}(count)=TempSpikeTime(t);        % ...keep the point - index it according to how many points have already been retained
           count=count+1; 
        elseif (TempSpikeTime(t)<(laststimpt-pre)-5) & (TempData(o,3+TempSpikeTime(t))<threshold) ; 
           TempSpikeTimeIndex1{o}(count)=TempSpikeTime(t);        
           count=count+1; 
        elseif (TempSpikeTime(t)<(laststimpt-pre)-5) & (TempData(o,4+TempSpikeTime(t))<threshold) ; 
           TempSpikeTimeIndex1{o}(count)=TempSpikeTime(t);        
           count=count+1; 
       elseif (TempSpikeTime(t)<(laststimpt-pre)-5) & (TempData(o,5+TempSpikeTime(t))<threshold) ; 
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
        [temp1(spike),temp2(spike)]=max(TempData(b,TempSpikeTimeIndex1{b}(1:transition(spike))));     % determine the time at which the spike reaches its negative peak
        TempSpikeTimeIndex2{b}(spike)=TempSpikeTimeIndex1{b}(temp2(spike))/10;                               % put the time, converted back to ms, in the index  
    end
    for spike=2:length(transition)                                                                       % for each putative spike except the last
        [temp1(spike),temp2(spike)]=max(TempData(b,TempSpikeTimeIndex1{b}(transition(spike-1)+1:transition(spike))));
        TempSpikeTimeIndex2{b}(spike)=TempSpikeTimeIndex1{b}(transition(spike-1)+temp2(spike))/10;
    end
    for spike=length(transition)+1
        [temp1(spike),temp2(spike)]=max(TempData(b,TempSpikeTimeIndex1{b}(transition(spike-1)+1:end)));
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


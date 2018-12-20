function [SpikeTrain,PSTH]=ConstructPSTH_JC(CellParameters,stimpoints) ;



% stimpoints= 95000;
tailpoints=0;

% construct a PSTH
clear SpikeTrain

for i=1:length(CellParameters);                      % for each of the trials at this mean light intensity....
    SpikeTrain(i,:)=zeros(1,(stimpoints+tailpoints));             % ...make a new matrix composed of zeros
    for j=1:length(CellParameters{i});                % for each of the spikes in this particular trial...
        a=CellParameters{i}(j);                    % ....at what time did the spike occur
        SpikeTrain(i,a*10)=1;                           % at the time - times 10, so that it can be index properly given the 0.1 ms sampling interval - the spike occurred replace the 0 with a 1
        clear a;                                                        
    end
end

[b,c]=size(SpikeTrain);                             % calculate the number of epochs represented in the SpikeTrain
PSTH=zeros(1,stimpoints+tailpoints);                   % create a new variable into which Spikes from different form different epochs will be binned

for k=1:b;                                      % for each of the epochs....
    PSTH=PSTH+SpikeTrain(k,:);                  % ....place into PSTH a value of 1 at the time at which a spike occurred 
end


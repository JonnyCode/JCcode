function [spks]=spks2(SpikeTimeIndex)

% JC 1/12/07
% this function takes Spike times from SpiketimeIndex and puts them into a
% matrix with times fo they can be ploted as a function of time. 


time=[0:.1:500]';
for d=1:length(SpikeTimeIndex);                     % for d from 1 to the number of cells in SpikeTimeIndex...
    spks{d}=zeros(1,5001)';                         % create a zero matrix for each cell in SpikeTimeIndex
    for b=1:length(SpikeTimeIndex{d});              % for b from 1 to the length of SpikeTimeIndex cell d... 
        c=(SpikeTimeIndex{d}(1,b).*10)+1.;          % create a matrix c equal to the b column of SpikeTimeIndex (cell d) times 10 
        spks{d}(c,1)=1;                             % make spks (cell d) column row c equal to 1
    end
    spks{d}=[time,spks{d}];
end
end
    


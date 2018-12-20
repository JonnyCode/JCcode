function [spkdpairs_norm,spkdpairs,spkdmean]=pairmaker_norm_FMR(SpikeTimeIndex)

%   This function makes all possible pairs out of the ephochs in 
% SpikeTimeIndex, runs the spkd function (distance metric) on each pair and 
% calculates the mean of of the spkd values.  

%JC 1/15/07

%   Additionally this function takes calculates normalized spike distance 
% for each cost from 0:1 in .01 cost intervals

% JC 2/6/07

% Create a spkdpairs{cell} for each cost which has all pairwise spkd comparisons

%close(figure_201, figure_200)

cost = [1:1:20];
for c=1:length(cost);                                                         % for each cost point in the spike distance metric     
    spkdpairs{c}=NaN((length(SpikeTimeIndex)-1),length(SpikeTimeIndex));   % Creates a matrix of NaNs
    for a=1:(length(SpikeTimeIndex)-1);                                    % each epoch a compared with epoch b 
        for b=a+1:length(SpikeTimeIndex);                                  % a and b (line below) allow for comparisons that have not been made(eg. 1v2, 1v3, 1v4, 2v3,2v4,3v4)
            d=spkd_c(SpikeTimeIndex{a},SpikeTimeIndex{b},length(SpikeTimeIndex{a}), length(SpikeTimeIndex{b}), cost(c));              % runs spike distance function at a given cost 
%            d=spkd(SpikeTimeIndex{a},SpikeTimeIndex{b}, cost);              % runs spike distance function at a given cost 
            spkdpairs{c}(a,b)=d ;                                           % sets the comparison of epoch a vs b into row a column b of matrix spkdpairs
        end
    end
end

%   Normalize across all cells so the highest spike distance for a given pair 
% across costs equals 1 and the lowest equals zero

for f=1:length(spkdpairs) ;                                                % for each cost cell
    spkdpairs_norm{f}= spkdpairs{1}*0 ;                                    % makes a matrix of approprate size
    for g=1:length(SpikeTimeIndex)-1 ;                                     % for each epoch g  
        for h=g+1:length(SpikeTimeIndex) ;                                 % compared with epoch h
        spkdpairs_norm{f}(g,h)=(spkdpairs{f}(g,h)-spkdpairs{1}(g,h))/(spkdpairs{length(spkdpairs)}(g,h)-spkdpairs{1}(g,h)) ;      % normalize each point in each cost matrix by subtracting the difference in spike number when cost=0 and dividing by the newly subtracted peak spike distance
        end
    end
end

% plots the first non-self-comparison across costs
for i=1:length(spkdpairs_norm) ;                                            % for each cost matrix
    plt(1,i)=spkdpairs_norm{i}(1,2) ;                                       % make a vector of spk distance values of epoch 1 vs epoch 2 
end
    figure(200) ;
    title('epoch 1 vs 2')
    xlabel('cost')
    ylabel('spike distance')
    axis([0 21 0 1.5]) ;
    plot(cost, plt,'r-') ;                                                  % plot distance value vs cost


% calculates the average of all distance metrics
for j=1:length(spkdpairs_norm)                                              %for each cost 
    spkdmean{j}=nanmean(spkdpairs{j}(1:end));                               % calcuate the spike distance mean        
end

% normalizes spkdmean

for k=1:length(spkdmean) ;                                            % for each cost matrix
    plt2(1,k)=spkdmean{k}(1,1) ;                                       % make a vector of average spk distance values  
end

spkdmean_norm=(plt2-min(plt2))/(max(plt2)-min(plt2))                    % normalizes to 1

% plot average spike distance against cost

figure(201) ;
axis([0 21 0 1.5]) ;
plot(cost, spkdmean_norm,'r-') ; 
title('average spike distance')
xlabel('cost')
ylabel('spike distance')




end


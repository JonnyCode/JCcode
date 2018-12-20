function [spkdpairs_norm,spkdmean,spkdmean_norm]=pairmaker_norm(spkdpairs)



%   Normalize across all cells so the highest spike distance for a given pair 
% across costs equals 1 and the lowest equals zero

cost=[1:1:20]

for f=1:length(spkdpairs) ;                                                % for each cost cell
    sub1 = (spkdpairs{f}-spkdpairs{1})                                      % subtract the min spike distance 
    sub2 = (spkdpairs{length(spkdpairs)}-spkdpairs{1}) ;                    % subtract the min spike distance from the max
    spkdpairs_norm{f}=sub1./sub2                                           % divide these two
end

% plots the first non-self-comparison across costs
for i=1:length(spkdpairs_norm) ;                                            % for each cost matrix
    plt(1,i)=spkdpairs_norm{i}(1,2) ;                                       % make a vector of spk distance values of epoch 1 vs epoch 2 
end

figure (202);
axis([0 21 0 1.5]) ;
plot(cost, plt,'g-') ;                                                  % plot distance value vs cost
title('epoch 1 vs 2')
xlabel('cost')
ylabel('spike distance')
hold on

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

figure (203);
axis([0 21 0 1.5]) ;
plot(cost, spkdmean_norm,'g-') ; 
title('average spike distance')
xlabel('cost')
ylabel('spike distance')
hold on



end


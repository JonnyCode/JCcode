function [spkdpairs]=pairmaker(SpikeTimeIndex)

%   This function makes all possible pairs out of the ephochs in 
% SpikeTimeIndex, runs the spkd function (distance metric) on each pair at
% a given cost.  The result is a matrix of spike distances at each given
% cost.


cost = [1:1:20];
for c=1:length(cost);                                                         % for each cost point in the spike distance metric     
    spkdpairs{c}=NaN((length(SpikeTimeIndex)-1),length(SpikeTimeIndex));      % Creates a matrix of NaNs
    for a=1:(length(SpikeTimeIndex)-1);                                       % each epoch a compared with epoch b 
        for b=a+1:length(SpikeTimeIndex);                                     % a and b (line below) allow for comparisons that have not been made(eg. 1v2, 1v3, 1v4, 2v3,2v4,3v4)
            d=spkd_c(SpikeTimeIndex{a},SpikeTimeIndex{b},length(SpikeTimeIndex{a}), length(SpikeTimeIndex{b}), cost(c));              % runs spike distance function at a given cost 
            spkdpairs{c}(a,b)=d ;                                            % sets the comparison of epoch a vs b into row a column b of matrix spkdpairs
        end
    end
end

end
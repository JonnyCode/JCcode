function rankedSignal = ranker(signal) ; 

% this function takes a matrix with each signal in rows and ranks the elements according to
% magnitude and in the case of a tie it makes the rank of the elements
% equal to the average of the ranks had those elements been slightly
% different, so that the sum of the ranked vector is equal to the
% sum([1:length(vector)]).  NaNs remain where they were.  
% This makes Nonparametric ranked correlation possible.

% JC 1/8/07

% signal = [-1,2,2,1,2,2,2,5.5,6 ; 1,2,1,2,2,5,6,7,7] ; % example signal to be ranked

[sorted, indexsorted] = sort(signal,2) ; % sort the signal into order of magnitude (sorted) and index of order of magnitude (indexsorted) 
ranks = repmat([1:length(sorted)],size(signal,1),1) ; % the ranks of each value in a matrix 

for b = 1:size(signal,1) ; % for each row in the matrix 
    numNaNs = sum(isnan(sorted(b,:)),2) ; % how many NaNs are there in the signal (NaNs get sorted last)
    ranks(b,end-numNaNs+1:end) = NaN ; % make the rank of NaNs a NaN
    for a = 1:length(sorted)-1 % for each value of sorted besides the last value ... 
        if sorted(b,a)==sorted(b,a+1) % if the next value is the same...
        endString = find([diff(sorted(b,a:end)),1]~=0,1) + a - 1 ; % find the first indice past "a" that is not the same number (adding a and subtracting 1 is because we are only looking at indices past a) also, the comma 1 after diff allows it to detect the last pair as a tie
        ranks(b,a:endString) = sum(ranks(b,a:endString))/(endString+1-a) ; % find the average rank of the string of repeated numbers
        end
    end
    rankedSignal(b,indexsorted(b,:)) = ranks(b,:) ; % index the new ranks vector by the order they came from the orignal vector
end


  
  






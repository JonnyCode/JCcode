% understanding sequence of LN models
% question 1: can a single LN model reproduce the exact output an sequence
% of LN models?
% question 2: is it possible to extract the LN model cascade by observing
% only the stimulus and final output.

tp = 10000 ; % time points

% LN filters
LinF1 = ; % linear filter 1
NonLinF1 = ; % nonlinear filter 2
LinF2 = ; 
NonLinF2 = ; 

% stimulus
stim = normrnd(0,10,1,tp) ; % stimulus
stim = lowPassFilter(stim) ;

% LN models
rspL1 = conv(stim,LinF1) ;
rspLN1 = ;
rspL2 = ;
rspLN2 = ;

% LN between stim and final response
LinF = ;
rspL = ;

figure % answer to question 1
plot(rspLN2, rspL,'*')

%













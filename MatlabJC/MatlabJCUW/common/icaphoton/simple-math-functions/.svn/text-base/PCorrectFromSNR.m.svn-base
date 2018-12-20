% return PCorrect for given SNR assuming Gaussian
% distribution
function PCorrect = PCorrectFromSNR(SNR)

% use MATLABs erf routine
PCorrect = 0.5 + erf(SNR/sqrt(2)) / 2;

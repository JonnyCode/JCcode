function SNR = SNRFromPCorrect(PCorrect, MaxSNR)

SNR = sqrt(2) * erfinv((PCorrect - 0.5) * 2);

% check for PCorrect = 1 (leads to infinite SNR)
for pnt = 1:length(PCorrect)
	if (PCorrect(pnt) >= 1)
		SNR(pnt) = MaxSNR;
	end
end

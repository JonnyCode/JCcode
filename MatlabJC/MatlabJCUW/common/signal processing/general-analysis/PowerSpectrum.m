function PowerSpec = PowerSpectrum(EpochData, PrePts, PSPts, SamplingInterval)
% PowerSpec(EpochData, PrePts, PSPts, SamplingInterval)
% 
% Compute power spectrum of data in EpochData matrix.  Start at PrePts and end
% at PrePts + PSPts.  SamplingInterval used to scale to power spectral
% density (i.e. power per unit bandwidth).
%
%%%SU
% load PowerSpectrumTest/test
% PS = PowerSpectrum(EpochData, PrePts, PSPts, SamplingInterval);
%
%%%TS isequal(PS, TargetPS)

[NumIterations, EpochPts] = size(EpochData);

for cnt = 1:NumIterations
	clear temp;
	temp(1:PSPts) = EpochData(cnt,(PrePts+1):(PrePts + PSPts));
	tempfft = fft(temp);
	if (cnt == 1)
		PowerSpec = tempfft .* conj(tempfft);
	else
		PowerSpec = PowerSpec + tempfft .* conj(tempfft);
	end
end

PowerSpec = real(PowerSpec) / NumIterations;
PowerSpec = 2 * PowerSpec * SamplingInterval / PSPts;

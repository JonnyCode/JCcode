% Linear approximation to differential equations describing 
% transduction cascade in vertebrate rods.  Standard
% set of differential equations describing cascade
% are linearized and solved.  Program computes
% linear filter approximation to transduction cascade and 
% convolves filter with rhodopsin time course vector provided
% as input to generate modeled response. (see e.g. Rieke and Baylor,
% 1998 RMP for differential equations)
%
% Created 7/01 FMR
% Revised 10/9/05 FMR
%   fixed bug that created a little structure at end of filter as negative
%   frequencies were not handled correctly.


% NOTE: what happened to 0 frequency term?

function [ModelResponse, Filter] = LinearCascadeModel(CascadeParameters, RhTimeCourse)

% dark cGMP concentration (in uM)
cGMPDark = ((CascadeParameters.DarkCurrent / CascadeParameters.cGMP2Cur).^(1/3));

% useful constant
con1 = 3 * CascadeParameters.CoopGC * CascadeParameters.PDEDark * CascadeParameters.ExRate;

% frequency vector in radians/sec
freq = 1:length(RhTimeCourse);
FreqStep = 1/(CascadeParameters.TimeStep * length(RhTimeCourse));
freq = 2 * pi * (freq-1) * FreqStep;

% generate linear filter approximation to impulse response of cascade
%from eq. 25, Rieke, Baylor 1996
filt = -cGMPDark ./ (CascadeParameters.PDEDecay + i .* freq);
temp = CascadeParameters.PDEDark + con1 * CascadeParameters.ExRate ./ (CascadeParameters.ExRate^2 + freq.^2);
temp = temp + i .* freq .* (1 - con1 ./ (CascadeParameters.ExRate^2 + freq.^2));
filt = filt ./ temp;

% get negative frequencies
for pnt = 1:length(RhTimeCourse)/2-1
	filt(length(RhTimeCourse) - pnt + 1) = conj(filt(pnt+1));
end

% calculate predicted flash response
rh = fft(RhTimeCourse);
cgmp = rh .* filt;
Filter = real(ifft(filt));
cgmp = real(ifft(cgmp));
ModelResponse = 3 * CascadeParameters.cGMP2Cur * cGMPDark^2 * cgmp;

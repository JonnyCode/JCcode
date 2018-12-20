function Filter = simFilter(time,tpeak,peakRise,peakAmp,ttrough,troughDecay,troughAmp)

% this function will approximate a linear filter by summing two gaussians
% and mulitplying negative gaussian by saturating exponetial to prevent negative values before peak

% JC 4/17/09

% time=.001:.001:.3 ;    % time vector
%
% tpeak = .05 ;   % time of peak
% peakRise = .010 ;   % rise of peak
%
% ttrough = .11 ;   % time of trough
% troughDecay = .06 ;    % decay of trough
%
% peakAmp = 1 ;    % sign and amplitude of peak
% troughAmp = .2 ;   % sign and amplitude of trough

Tr = tpeak ;    % 

y1 = exp(-((time-tpeak).^2)/(2*peakRise^2)) ;  % make gaussian for peak
y2 = exp(-((time-ttrough).^2)/(2*troughDecay^2)) ;  % make gausian for trough
y3 = exp((-Tr./time).^3) ;                 % make sat exponential to avoid tough before peak

Filter = (peakAmp*y1)-((troughAmp*y2).*y3) ;    % create filter

% figure
% plot(time,peakAmp*y1)
% hold on
% plot(time,-troughAmp*y2)
% plot(time,-troughAmp*y3)
% plot(time,Filter,'r--')

end


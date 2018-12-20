function omega = FreqAxis(N)
% omega = FreqAxis(EpochCondition(cond))
%
% function returns a vector of frequencies in Hz
% This is for plotting symetrical transform plots so you want symetrical
% frequencies as well.
% len specifies the length of the epoch.

% Number of steps gives number of ms in the trial: 1000 steps give 1 sec
deltaOmega = 1/N*1000;
omega1 = [0:deltaOmega:deltaOmega*N/2];
omega2 = -fliplr(omega1(2:end));
omega = [omega2,omega1(1:end-1)];
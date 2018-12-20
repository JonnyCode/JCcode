function t = TimeAxis(N)
%  t = TimeAxis(N)
%
% function returns a vector of times in sec
% N specifies the number of time points of the epoch.

% Number of steps gives number of ms in the trial: 1000 steps give 1 sec

deltaT = 1/1000;
t = [0:deltaT:N/1000-deltaT];
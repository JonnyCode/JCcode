function[null] = nullmin(rho)

%Function returns spike number of direction with least firing

% Sneha Ravi 
% Last revision: 12-18-2012

null = cell(length(rho), 1);%
for i = 1:length(rho)
    null{i,1} = min(rho{i,1},[], 2)';
end
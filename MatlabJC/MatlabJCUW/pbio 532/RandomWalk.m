function walk = RandomWalk(StepSize, NumSteps)

% Core function for random walk simulations.  Create walk for specified
% number of steps and step size, assuming steps independent, equally likely
% positive or negative, and equal length. 

walk(1) = 0;
walk(2:NumSteps+1) = cumsum((floor(unidrnd(2, NumSteps, 1)/2)*2-1) * StepSize);


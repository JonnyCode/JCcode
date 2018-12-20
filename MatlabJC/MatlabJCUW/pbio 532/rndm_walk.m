%--------------------------------------------------------------
% RandomWalkScript
%   Demonstrate several important properties of random walks in biology.
%
% Contents:
%   First part of code (up to line 150 or so) contains simulation I used
%   for part of lecture in class.  Second part (beyond line 150) contains
%   code and questions for class discussion. 
%
% Created 9/03 FMR for NeuBeh/PBio 532
% Last Revised 12/06 FMR
%   - general cleaning
%   - added discussion problems
%--------------------------------------------------------------

%--------------------------------------------------------------
% Simple examples of random walks
%   RandomWalk function simulates steps to right or left of equal probability.
%   Steps are equal in size.  
%--------------------------------------------------------------

% plot random walk one step at a time
StepSize = 1;
clear walk;
figure(1);

% repeat until told to halt (0 entered)
while (1)
    if (exist('walk'))
        walk = [walk, RandomWalk(StepSize, 1)];     % append one step to existing trajectory
    else
        walk = RandomWalk(StepSize, 1);             % start new trajectory
    end
    clf;
    axes('FontSize', 18);
    plot([1:length(walk)], cumsum(walk));           % random trajectory
    hold on;
    plot([1 length(walk)], [0 0], 'r');             % horizonal line at start position
    hold off;
    ylabel('x-position');
    xlabel('time');
    title('1-D random walk', 'FontSize', 18);
    axis([0 length(walk) min(cumsum(walk))-0.2 max(cumsum(walk))+0.2 ])
    KeyEvent = input('0 to quit\n');                % continue or halt?
    if (KeyEvent == 0)
        break;
    end
end    

% plot example random walks with specified number of steps
figure(1);
steps = input('Input number of steps\n');
clear walk;
while (1)
    walk = RandomWalk(StepSize, steps);
    clf;
    axes('FontSize', 18);
    plot([1:length(walk)], walk);
    hold on
    plot([1 length(walk)], [0 0], 'r');
    hold off
    ylabel('x-position');
    xlabel('time');
    title('1-D random walk','FontSize', 18);
    axis([0 length(walk) min(walk)-0.2 max(walk)+0.2 ])
    KeyEvent = input('0 to quit\n');
    if (KeyEvent == 0)
        break;
    end
end

%--------------------------------------------------------------
% Statistical properties
%--------------------------------------------------------------

% distribution of final positions for walk of specified length
figure(1); clf;
axes('FontSize', 18);
Trials = 4000;          % number of walks
steps = input('Input number of steps\n');
clear walk;
clear FinalPos;
for trial = 1:Trials
    walk = RandomWalk(StepSize, steps);
    FinalPos(trial) = walk(steps+1);
end

% rms position should be proportional to the sqrt of the number of steps
fprintf(1, 'RMS position = %d\n', std(FinalPos));
hist(FinalPos, 20);
xlabel('final position');
ylabel('count');
title('Distribution of final positions', 'FontSize', 18);
xlim([-100 100])

%--------------------------------------------------------------
% Diffusion in confined spaces (e.g. inside cellular compartment)
%--------------------------------------------------------------

% parameters
DiffusionCoefficient = 1;           % um^2/msec (typical ion)
TimeStep = 0.1;                     % in msec
StepSize = sqrt(2*DiffusionCoefficient*TimeStep);     % in um
CellSizeX = 10;                   % x cell size in um
CellSizeY = 10;                   % y cell size in um

% initial conditions
clear walk;
steps = 5000;
walkx(1) = 0;
walky(1) = 0;

figure(1);
clf;
axes('FontSize', 18);

% unconfined random walk in x-y space, x and y independent
walkx(2:steps+1) = cumsum(normrnd(0, StepSize, 1, steps));
walky(2:steps+1) = cumsum(normrnd(0, StepSize, 1, steps));
for step = 1:steps
    % check whether walk exceeds boundary in each direction, if so assume
    % reflection from wall
    if (walkx(step) > CellSizeX/2)
        walkx(step:steps) = walkx(step:steps) - 2 * (walkx(step) - CellSizeX/2);
    end
    if (walkx(step) < -CellSizeX/2)
        walkx(step:steps) = walkx(step:steps) - 2 * (walkx(step) + CellSizeX/2);
    end
    if (walky(step) > CellSizeY/2)
        walky(step:steps) = walky(step:steps) - 2 * (walky(step) - CellSizeY/2);
    end
    if (walky(step) < -CellSizeY/2)
        walky(step:steps) = walky(step:steps) - 2 * (walky(step) + CellSizeY/2);
    end
    % update plot of trajectory every 50 steps
    if (rem(step, steps/50) == 0)
        plot(walkx(1:step), walky(1:step));
       	hold on;
        plot([-CellSizeX/2 CellSizeX/2 CellSizeX/2 -CellSizeX/2 -CellSizeX/2], [-CellSizeY/2 -CellSizeY/2 CellSizeY/2 CellSizeY/2 -CellSizeY/2], 'r');
        AxisLimits = max([CellSizeX CellSizeY]);
        axis([-AxisLimits/2*1.1 AxisLimits/2*1.1 -AxisLimits/2*1.1 AxisLimits/2*1.1]);
        hold off
		xlabel('x-position (um)');
		ylabel('y-position (um)');
		title('2-D random walk', 'FontSize', 18);
        TextString = strcat(num2str(step * TimeStep),' msec');
        text(CellSizeX*0.3, CellSizeX*0.4, TextString, 'FontSize', 16);
        pause(0.2)
    end
end    

%--------------------------------------------------------------
% Problem Set #1: Diffusion models and spiking statistics
%   One common model for a spiking cell is that the subthreshold membrane
%   potential behaves as a random walk due to a barrage of excitatory and
%   inhibitory inputs, producing positive and negative steps in membrane
%   potential.  One nice aspect of such models is that the analytical tools
%   for understanding random walks applies to the dynamics of the
%   subthreshold voltages.  Another is that such models are easy and fast
%   to simulate.  The code below explores a simple version of such a model.
%
%   This code simulates a cell that receives random excitatory and
%   inhibitory inputs.  Excitatory inputs cause depolarizing steps and
%   inhibitory inputs hyperpolarizing steps.  When the resulting random
%   walk reaches a threshold an action potential is triggered and the
%   position of the walk reset to 0.
%--------------------------------------------------------------

% parameters
ExcitatoryProb = 0.5;               % probability of making excitatory step
InhibitoryProb = 1-ExcitatoryProb;  % probability of making inhibitory step
ExcitatoryStepSize = 1.5;             % size of excitatory step
InhibitoryStepSize = -1;            % size of inhibitory step
Threshold = 10;                     % spike threshold
steps = 20000;                     % number of time steps to simulate

% random walk without thresholding, independent excitatory and inhibitory
% steps with appropriate step sizes
clear walk SpikeTrain;
walk(1) = 0;

% first random walk without worrying about step sizes
walk(2:steps) = unifrnd(-InhibitoryProb, ExcitatoryProb, steps-1, 1);
% now set sizes of negative and positive steps
Indices = find(walk < 0);
walk(Indices) = InhibitoryStepSize;
Indices = find(walk > 0);
walk(Indices) = ExcitatoryStepSize;
walk = cumsum(walk);

% apply threshold and resets
SpikeTrain(1:steps) = 0;
for pnt = 1:steps
    if (walk(pnt) > Threshold)
        SpikeTrain(pnt) = 1;
        walk(pnt+1:steps) = walk(pnt+1:steps) - walk(pnt);
    end
end

% plot random walk and spike times
figure(1); clf;
plot(walk);
hold on;
plot(SpikeTrain * std(walk) + Threshold, 'r');
hold off;
xlim([0 1000]);
xlabel('time');
ylabel('position');

% inter-spike interval statistics
SpikeTime = find(SpikeTrain == 1);
Intervals = SpikeTime(2:length(SpikeTime)) - SpikeTime(1:length(SpikeTime)-1);
figure(2);
hist(Intervals, 50);
xlabel('interval');
ylabel('count');

% (1) Why does the inter-spike interval histogram have the shape that it
% does?  Why are there no short intervals - i.e. what are the analogs of
% the absolute and relative refractory periods?  Can you test your proposed
% explanations?
%
% (2) How does the shape of the interval distribution depend on the step
% sizes and probabilities for excitatory and inhibitory inputs? Can you
% explain the changes intuitively?
%
% (3) What happens if the excitatory and inhibitory steps are equal in size
% and probability?  Is this realistic biologically?  What options (in
% general - not just in code above) do you have to fix it?

%--------------------------------------------------------------
% TO HERE FOR 1/22 DISCUSSION
%--------------------------------------------------------------

%--------------------------------------------------------------
% Comparison of random walk predictions with Fick's equations for
% macroscopic diffusion
%--------------------------------------------------------------

% distribution of final positions for random walk compared to solution of
% Fick's equations for point source

DiffusionCoefficient = 0.1;           % um^2/msec (typical ion)
TimeStep = 0.001;                     % in msec
StepSize = sqrt(2*DiffusionCoefficient*TimeStep);     % in um
Trials = 1000;
steps = 5000;

figure(1); clf;
axes('FontSize', 18);
clear walk FinalPos;
for trial = 1:Trials
    walk = RandomWalk(StepSize, steps);
    FinalPos(trial) = walk(steps+1);
end
[dist, distx] = hist(FinalPos, 20);
dist = dist / Trials;
bar(distx, dist);
xlabel('final position (um)');
ylabel('probability');
title('Distribution of final positions', 'FontSize', 18);
xlim([-std(FinalPos)*4 std(FinalPos)*4])
hold on

% predicted distribution from Fick's Equations
predict = (distx(2) - distx(1)) * exp(-distx.^2/(4*DiffusionCoefficient*steps * TimeStep)) / sqrt(4 * 3.14159 * DiffusionCoefficient * steps * TimeStep);
plot(distx, predict, 'LineWidth', 2, 'Color', [1 0 0]);
legend('random walk','predicted from Ficks eqns');
legend boxoff;

% rms position
fprintf(1, 'RMS position from random walk = %d\n', std(FinalPos));
fprintf(1, 'RMS position from macroscopic diffusion = %d\n', sqrt(2 * DiffusionCoefficient * steps * TimeStep));

%--------------------------------------------------------------
% Diffusion-limited reactions (in 2D)
%   The rate of diffusional encounters between molecules sets a fundamental
%   limit to the speed of biochemical interactions.  The point of this
%   section of code is to explore how this limit depends on geometry,
%   concentration, size, etc.  
%
%   A specific example here is phototransduction: the initial events in the
%   conversion of light into an electrical signal take place on membrane
%   discs within the photoreceptor.  The speed of this phototransduction
%   process cannot exceed that set by the rate at which the light-activated
%   photopigment molecule encounters G-proteins.  The geometry and other
%   parameters below are a rough approximation of the situation on a single
%   disc.
%--------------------------------------------------------------

% parameters
DiffusionCoefficient = 0.1;           % um^2/msec (typical ion)
TimeStep = 0.001;                     % in msec
StepSize = sqrt(2*DiffusionCoefficient*TimeStep);     % in um
CellSizeX = 1.5;                    % x cell size in um
CellSizeY = 1.5;                    % y cell size in um
NumEffectors = 1000;                % total number of effectors in cell
InteractionRadius = 0.01;               % in um

% initial conditions
clear walk;
steps = 400;
walkx(1) = 0;
walky(1) = 0;
clear Effectorx Effectory walkx walky
Effectorx(1:NumEffectors) = unifrnd(-CellSizeX/2, CellSizeX/2, NumEffectors, 1);
Effectory(1:NumEffectors) = unifrnd(-CellSizeY/2, CellSizeY/2, NumEffectors, 1);

% unconfined random walk in x-y space, x and y independent - Gaussian
% because each step here represents many discrete steps in microscopic
% random walk
walkx(2:steps+1) = cumsum(normrnd(0, StepSize, 1, steps));
walky(2:steps+1) = cumsum(normrnd(0, StepSize, 1, steps));

NumberActivations = 0;
AlreadyActivated = NaN;
for step = 1:steps
    % check whether walk exceeds boundary in each direction, if so assume
    % reflection from wall
    if (walkx(step) > CellSizeX/2)
        walkx(step:steps) = walkx(step:steps) - 2 * (walkx(step) - CellSizeX/2);
    end
    if (walkx(step) < -CellSizeX/2)
        walkx(step:steps) = walkx(step:steps) - 2 * (walkx(step) + CellSizeX/2);
    end
    if (walky(step) > CellSizeY/2)
        walky(step:steps) = walky(step:steps) - 2 * (walky(step) - CellSizeY/2);
    end
    if (walky(step) < -CellSizeY/2)
        walky(step:steps) = walky(step:steps) - 2 * (walky(step) + CellSizeY/2);
    end
    % check whether x,y position matches any of effector
    Distances = sqrt((walkx(step) - Effectorx).^2 + (walky(step) - Effectory).^2);
    Indices = find(Distances < InteractionRadius);
    if (~isempty(Indices))
        for effector = 1:length(Indices)
            SecondIndices = find(Indices(effector) == AlreadyActivated);
            if (isempty(SecondIndices))
                NumberActivations = NumberActivations + 1;
                AlreadyActivated(NumberActivations) = Indices(effector);
            end
        end
    end
end    

figure(1); clf;
axes('FontSize', 18);
plot(walkx, walky);
hold on;
plot([-CellSizeX/2 CellSizeX/2 CellSizeX/2 -CellSizeX/2 -CellSizeX/2], [-CellSizeY/2 -CellSizeY/2 CellSizeY/2 CellSizeY/2 -CellSizeY/2], 'r');
AxisLimits = max([CellSizeX CellSizeY]);
axis([-AxisLimits/2*1.1 AxisLimits/2*1.1 -AxisLimits/2*1.1 AxisLimits/2*1.1]);
xlabel('x-position (um)');
ylabel('y-position (um)');
title('2-D random walk', 'FontSize', 18);
plot(Effectorx, Effectory, '.k');
plot(Effectorx(AlreadyActivated), Effectory(AlreadyActivated), '.r');
hold off

% activation rate in 1/sec
fprintf(1, 'number activated = %d activation rate = %d\n', NumberActivations, NumberActivations / (steps * TimeStep * 0.001));

% (1) Why can we use a gaussian random variable for the random walk above,
% when our microscopic random walk picture is for discrete steps?  Why is
% the step size sqrt(2*DiffusionCoefficient*TimeStep)?
%
% (2) Some biochemical reactions take place on a surface - the cell surface
% membrane, or in the case of phototransduction an internal membrane.  Does
% the diffusion limit depend on the number of dimensions (assuming all else
% equal - e.g. diffusion coefficient, interaction radius, number of
% effectors, ...)?  Generalize the simulation above to three dimensions
% (assume a cubical cell) and determine the relative rates in 2D and 3D.
%
% (3) How does the activation rate depend on diffusion coefficient,
% interaction radius, concentration of effectors?  Why?
%
% (4) What could cause a given reaction to fall short of the diffusion limit?
%
% (5) What does the interaction radius represent physically?
%
% (6) Play around with the conditions above and report anything interesting
% or unexpected that you find.



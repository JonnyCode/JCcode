%****************************************************************************
% Differential equations tutorial for NeuBeh/PBIO 545, Winter 2003
%
% Created 12/02 Fred Rieke
% Revisions: 1/05 FMR updated for 2005 class
%
%****************************************************************************

% useful default stuff
set(0, 'DefaultAxesFontName','Palatino')
set(0, 'DefaultAxesFontSize', 12)
colormap([0 0 0])

%****************************************************************************
% Example 1: second-messenger cascades and phototransduction
%****************************************************************************
% The examples below build up to a model for phototransduction based
% on known chemical interactions.  We start with some of the building
% blocks for a more complicated model and work up to a full model for
% phototransduction.  Many of the examples below are much more generally
% applicable - e.g. the simple feedback models.  Phototransduction is a 
% nice example because we can construct a model based on known mechanisms
% with measured properties.
%****************************************************************************

%----------------------------------------------------------------------------
% Consider creation of a substance x by another y.  For example, y could be 
% an active receptor and x its downstream effector.  The rate of creation of x 
% is proportional to the amount of active y and x decays with a rate 
% constant alpha.  The differential equation describing this situation is:
% 	dx/dt = y - alpha * x
% This equation does not have a unique solution unless we also specify an initial 
% condition.  This is an important point that applies to both attempts to 
% solve differential equations numerically or analytically. We'll use x=0 at time 0.

% The approach we will take to solving this equation is to turn the differential
% equation into a difference equation, which we'll solve for finite time steps.  In
% doing this we will assume that the behavior of x is determined entirely by 
% its first derivative for our finite time step.  This amounts to a Taylor series
% approximation to x(t), where we ignore the second derivative and higher terms.  
% The accuracy of the approximation improves as the time step becomes smaller.  
% If this is unclear come back to it after going through the program below.  

figure(1);
x(1) = 0;				% initial condition
alpha = 20;				% rate constant for decay of x (in 1/sec)
TimeStep = 0.001;		% time step for difference equation (in sec)
PrePts = 200;			% points before step in y
StmPts = 400;			% number of time points that y is active
NumPts = 1000;			% total points to simulate

% initialize y; in this case y is a simple step
y(1:PrePts) = 0;
y(PrePts+1:PrePts + StmPts) = 1;
y(PrePts + StmPts + 1:NumPts) = 0;

% plot y
tme = 1:NumPts;
tme = (tme - PrePts) * TimeStep;
subplot(2, 1, 1);
plot(tme, y);
ylim([-.1 1.1]);
xlabel('time (sec)');
ylabel('input y');
title('activation of x by y, no feedback');

% Calculate x by making the differential equation above into a difference equation.  
% Rather than solving with time as a continuous variable, we discretize time into 
% finite steps of length TimeStep.  We will approximate x in each time step.  Thus
% the derivative
%	 dx/dt
% becomes
% 	 [x(n) - x(n-1)] / TimeStep
% where x(n) is the value of x in the nth time bin (i.e. at time n*TimeStep).
% Note that in the limit where TimeStep goes to 0 this is the definition of
% a derivative.  Make sure you see this connection as it is at the core of 
% how differential equations are solved numerically.
% So now our original differential equation becomes
%	[x(n) - x(n-1)] / TimeStep = y(n-1) - alpha * x(n-1)
% or, solving for x(n),
%	x(n) = x(n-1) + TimeStep * [y(n-1) - alpha * x(n-1)]
% This gives us an update rule for x - given x and y in time bin n-1, we 
% can compute x in time bin n.  Iterating this update rule gives us 
% an approximation to the full time course of x:
for pnt = 2:NumPts
	x(pnt) = x(pnt-1) + (y(pnt-1) - alpha * x(pnt-1)) * TimeStep;
end
subplot(2, 1, 2);
plot(tme, x);
xlabel('time (sec)');
ylabel('output x');


% This is an important example because it is one you use all the time.  The
% analytical solution (see lecture notes) for this equation is an
% exponential.  There are lots of processes you will approximate as single
% exponentials (e.g. the membrane time constant of a cell).  

% Questions/Problems #1:
% (1) play with the rate constant alpha and explain what happens.
% (2) try different shaped inputs (other than step) and explain the results.
% (3) play with TimeStep to see over what range of time bins the numerical solution
%	  is accurate (i.e. compare to an analytical solution to the same DEQ).

%----------------------------------------------------------------------------
% Let's elaborate this example a little.  What if x also has a rate of spontaneous activation,
% in addition to activation by y?  Now our differential equation becomes
%	dx/dt = y + eta - alpha * x
% where eta is the rate of spontaneous activation.
% We'll solve this in the same way.

figure(2);
eta = 1;			% rate of spontaneous activation of x
x(1) = eta / alpha;	% new initial condition

% plot input y
subplot(2, 1, 1);
plot(tme, y);
ylim([-.1 1.1]);
xlabel('time (sec)');
ylabel('input');
title('activation of x by y, no feedback');

% calculate x by making differential equation into difference equation
for pnt = 2:NumPts
	x(pnt) = x(pnt-1) + (y(pnt-1) + eta - alpha * x(pnt-1)) * TimeStep;
end
subplot(2, 1, 2);
plot(tme, x);
xlabel('time (sec)');
ylabel('output');

% Make sure you understand the form of the difference equation in the loop above.

% Questions/Problems #2:
% (1) How does the solution differ in this case from that above?  Why?  
% (2) Can you explain the change from the differential equation?
% (3) Why is eta/alpha is reasonable initial condition?  

%----------------------------------------------------------------------------
% The examples above can both be solved analytically, providing a useful check to 
% the numerical solution.  Now let's add a feedback term, which makes the analytical
% solutions difficult.  So now active x will feedback to modify the 'effective' 
% activity of y (that is rate of creation of x by y).  Now the differential
% equation becomes
% 	dx/dt = y * (1+gx)^n - alpha * x,
% where g is a constant describing the gain of the feedback and and n is
% an exponent determining the linear or nonlinear behavior of the feedback signal.
% There are certainly other ways a feedback mechanism could work (and correspondingly
% different differential equations), but this is one reasonable form.

figure(3);
x(1) = 0;						% initial condition
Power = -1;						% feedback power
g = 20;							% feedback gain

% plot input y
subplot(2, 1, 1);
plot(tme, y);
ylim([-.1 1.1]);
xlabel('time (sec)');
ylabel('input');
title('activation of x by y, feedback');

% solve difference equation
for pnt = 2:NumPts
	x(pnt) = x(pnt-1) + (y(pnt) * (1 + g * x(pnt-1))^Power - alpha * x(pnt-1)) * TimeStep;
end
subplot(2, 1, 2);
plot(tme, x);
xlabel('time (sec)');
ylabel('output');

% Questions/Problems #3:
% (1) How does the behavior compare to the case without feedback?
%	  Compare both the amplitude and kinetics of x.
% (2) How does effect of feedback change as power changed?  Why? Try both negative
%	  and positive values for the power.
% (3) How does changing g change things?  Why?
% (4) What happens when you change TimeStep?  Is the numerical solution more 
%	  or less sensitive to changes in TimeStep than that above?  Why?

%----------------------------------------------------------------------------
% Now let's consider combinations of a couple steps.  Anticipating building up
% to a model for phototransduction we'll change from the generic variables x and y
% to r (rhodopsin activity) and p (phosphodiesterase activity).  Light produces
% an electrical signal in the photoreceptor by activating rhodopsin which then
% activates phosphodiesterase (through transducin).  Active phosphodiesterase 
% breaks down cGMP and reduces the membrane current through cGMP-gated channels.  
% Start by considering phosphodiesterase activation.  First, we'll assume that the 
% rhodopsin activity (its ability to activate phosphodiesterase through transducin) 
% obeys
%	dr/dt = -sigma*r
% where sigma is the rate constant for the decay of rhodopsin.  Then 
% we'll assume that the phosphodiesterase activity obeys
%	dp/dt = r + eta - phi * p 
% where eta represents spontaneous phosphodiesterase activation and
% phi is the rate constant for decay.

figure(1)
sigma = 5;				% rhodopsin activity decay rate constant (1/sec)
phi = 5;				% phosphodiesterase activity decay rate constant (1/sec)
eta = 10;				% phosphodiesterase activation rate constant (1/sec)
r(1) = 1;				% initial condition for r
p(1) = eta/phi;			% initial condition for p

NumPts = 1000;			% number of points to simulate
TimeStep = 0.001;		% time step 

tme = 1:NumPts;
tme = tme * TimeStep;

% solve difference equation
for pnt = 2:NumPts
	r(pnt) = r(pnt-1) + TimeStep * (-sigma * r(pnt-1));
	p(pnt) = p(pnt-1) + TimeStep * (r(pnt-1) + eta - phi * p(pnt-1));
end

% plot time course of rhodopsin activity
subplot(2, 1, 1);
plot(tme, r);
xlabel('time (sec)')
ylabel('rhodopsin activity')
% plot time course of phosphodiesterase activity
subplot(2, 1, 2);
plot(tme, p);
xlabel('time (sec)')
ylabel('pde activity')

% Questions/Problems #4:
% (1) why do we choose an exponential decay for the shape of rhodopsin's activity?  What might
%	  change that?
% (2) Explore different combinations of sigma and phi and their impact.

%----------------------------------------------------------------------------
% The code above describes the 'activation' part of phototransduction - i.e. how
% light activation of rhodopsin leads to activation of phosphodiesterase.  Now 
% lets add in the steps linking that to a change in current, and the steps 
% required for the light response to recover.  The role of activated phosphodiesterase
% is to hydrolyze cGMP.  Another enzyme, guanylate cyclase, synthesizes cGMP.  So the 
% cGMP concentration g is controlled by a balance of synthesis (at a rate s) and hydrolysis
% (at a rate pg):
%	dg/dt = s - pg
% The membrane current depends on the third power of the cGMP concentration
%	I = k g^3
% This is all we need for a simple phototransduction model.  Below we add a feedback
% that controls the rate of synthesis.

figure(1)
gdark = 15;				% concentration of cGMP in darkness
cgmp2cur = 8e-3;		% constant relating cGMP to current

tme = 1:NumPts;
tme = tme * TimeStep;
g(1) = gdark;							% initial condition
s(1:NumPts) = gdark * eta/phi;		% steady state: synthesis = hydrolysis
p(1) = eta/phi;

% solve difference equations for each component
for pnt = 2:NumPts
	r(pnt) = r(pnt-1) + TimeStep * (-sigma * r(pnt-1));
	p(pnt) = p(pnt-1) + TimeStep * (r(pnt-1) + eta - phi * p(pnt-1));
	g(pnt) = g(pnt-1) + TimeStep * (s(pnt-1) - p(pnt-1) * g(pnt-1));
end

% determine current change
cur = cgmp2cur * g.^3;

% plot current, pde activity, synthesis rate, cGMP concentration
subplot(4, 1, 1);
plot(tme, cur);
xlabel('time (sec)');
ylabel('current');
subplot(4, 1, 2);
plot(tme, p);
xlabel('time (sec)');
ylabel('pde activity');
subplot(4, 1, 3);
plot(tme, s);
xlabel('time (sec)');
ylabel('synthesis rate');
subplot(4, 1, 4);
plot(tme, g);
xlabel('time (sec)');
ylabel('[cGMP]');

%----------------------------------------------------------------------------
% Now let's add a calcium feedback to the model above.  Calcium enters the outer
% segment through the cGMP-gated channels.  Hence the rate of calcium influx
% is proportional to the current flowing.  Calcium is removed by a Na+/K+,Ca2+ 
% exchanger.  The rate of removal is proportional to the calcium concentration:
%	dc/dt = qI - beta c
% where c is the calcium concentration, q is the proportionality constant between
% changes in calcium and the current, and beta is the rate constant for calcium
% removal.  
% Calcium acts on the rate of cGMP synthesis:
%   s = smax / (1 + (c/K)^h)
% where smax is the maximum rate and K and H are constant (affinity and cooperativity)
% describing the strength of the feedback.  Otherwise the model is the same.

figure(2)
cdark = 0.5;			% dark calcium concentration
beta = 20;				% rate constant for calcium removal in 1/sec
hillcoef = 4;			% cooperativity
hillaffinity = 0.3;		% affinity

cur2ca = beta * cdark / (cgmp2cur * gdark^3);				% get q using steady state
smax = eta/phi * gdark * (1 + (cdark / hillaffinity)^hillcoef);		% get smax using steady state

tme = 1:NumPts;
tme = tme * TimeStep;

% initial conditions
g(1) = gdark;
s(1) = gdark * eta/phi;		
c(1) = cdark;
p(1) = eta/phi;

% solve difference equations
for pnt = 2:NumPts
	r(pnt) = r(pnt-1) + TimeStep * (-sigma * r(pnt-1));
	p(pnt) = p(pnt-1) + TimeStep * (r(pnt-1) + eta - phi * p(pnt-1));
	c(pnt) = c(pnt-1) + TimeStep * (cur2ca * cgmp2cur * g(pnt-1)^3 - beta * c(pnt-1));
	s(pnt) = smax / (1 + (c(pnt) / hillaffinity)^hillcoef);
	g(pnt) = g(pnt-1) + TimeStep * (s(pnt-1) - p(pnt-1) * g(pnt-1));
end
% determine current change
cur = cgmp2cur * g.^3;

% plot current, pde, synthesis, cGMP and calcium
subplot(5, 1, 1);
plot(tme, cur);
xlabel('time (sec)');
ylabel('current');
subplot(5, 1, 2);
plot(tme, p);
xlabel('time (sec)');
ylabel('pde activity');
subplot(5, 1, 3);
plot(tme, s);
xlabel('time (sec)');
ylabel('synthesis rate');
subplot(5, 1, 4);
plot(tme, g);
xlabel('time (sec)');
ylabel('[cGMP]');
subplot(5, 1, 5) 
plot(tme, c)
xlabel('time (sec)');
ylabel('[calcium]');

% Questions/Problems #5:
% (1) Explain the relation between the steady state conditions
%	  and the constants q and smax.
% (2) Play with the various parameters and see how the alter the calculated 
%	  light response.  Explain why things change they way they do?
% (3) The model will generate damped oscillations for some values of beta.  Why?

%----------------------------------------------------------------------------
% Example 2: two-state systems and channel gating particles
%----------------------------------------------------------------------------

% The Hodgkin-Huxley model relies on a series of gating particles, each 
% of which can be active or inactive.  Lets consider one of them - the m gate for
% sodium channel opening. The probability that m is in the active state (m=1) 
% is p, and the probability it is in the inactive state (m=0) is (1-p).  
% The time derivative of m is given by the likelihood that it becomes active 
% (alpha * (1-p)) minus the likelihood that it inactivates (beta p).  
%	dm/dt = alpha(1-p) - beta p

figure(1)
m(1) = 1;				% initial condition
alpha = 200;				% active->inactive rate constant
beta = 200;				% inactive->active rate constant
NumPts = 1000;			
TimeStep = 0.00001;

% solve difference equation
for pnt = 2:NumPts
	m(pnt) = m(pnt-1) + TimeStep * ((1 - m(pnt-1)) * alpha - m(pnt-1) * beta);
end

% plot results
tme = 1:NumPts;
tme = tme * TimeStep;
plot(tme, m);

% This is what we might expect if we look back the the differential
% equation for m above.  We can rearrange things a little and see that the
% solution is an exponential dependence on time.

%----------------------------------------------------------------------------
% Now let's make alpha and beta voltage dependent

figure(2);
PrePts = 200;			% points before step in y
StmPts = 400;			% number of time points that y is active
m(1) = 0;				% initial condition

% initialize voltage - step from -60 to -40 in for 200 points
v(1:PrePts) = -60;
v(PrePts+1:PrePts + StmPts) = -40;
v(PrePts + StmPts + 1:NumPts) = -60;

% compute alpha from HH parameters
alpha(1:NumPts) = -100 * (v + 30) / (exp(-(v + 30)/10) - 1);

% compute beta from HH parameters
beta(1:NumPts) = 4000 * exp(-(v + 55) / 18);

% solve difference equation
for pnt = 2:NumPts
	m(pnt) = m(pnt-1) + TimeStep * ((1 - m(pnt-1)) * alpha(pnt-1) - m(pnt-1) * beta(pnt-1));
end

% plot results
tme = 1:NumPts;
tme = tme * TimeStep;
plot(tme, m);

% Some extensions (not problems to be turned in):
% (1) add an inactivation gating particle
% (2) build this up to the HH sodium current model

%----------------------------------------------------------------------------
% Fourier transforms and differential equations
%----------------------------------------------------------------------------

% Fourier transforms can provide a nice approach to solve differential
% equations.  The key trick is that the derivatives vanish when we take the
% Fourier transform - i.e. if c(t) and C(w) are Fourier transform pairs, so
% are dc(t)/dt and -2iwC(w).  Make sure you understand why this is true
% from the definition of the Fourier transform.

% Let's use this to solve the first differential equation we considered
% above
%    	dx/dt = y - alpha * x
% If we take the Fourier transform of both sides of this equation, we get
%       -iwX(w) = Y(w) - alpha * X(w) 
% So
%       X(w) = Y(w) / (alpha + i w)
% Let's make sure this works.  First, solve this as we did initially by
% constructing a difference equation:

figure(1);
x(1) = 0;				% initial condition
alpha = 20;				% rate constant for decay of x (in 1/sec)
TimeStep = 0.001;		% time step for difference equation (in sec)
PrePts = 200;			% points before step in y
StmPts = 400;			% number of time points that y is active
NumPts = 1000;			% total points to simulate

% initialize y; in this case y is a simple step
y(1:PrePts) = 0;
y(PrePts+1:PrePts + StmPts) = 1;
y(PrePts + StmPts + 1:NumPts) = 0;

% plot y
tme = 1:NumPts;
tme = (tme - PrePts) * TimeStep;
subplot(3, 1, 1);
plot(tme, y);
ylim([-.1 1.1]);
xlabel('time (sec)');
ylabel('input y');

% solve difference equation for x
title('activation of x by y, no feedback');
for pnt = 2:NumPts
	x(pnt) = x(pnt-1) + (y(pnt-1) - alpha * x(pnt-1)) * TimeStep;
end
subplot(3, 1, 2);
ylim([0 0.06]);
plot(tme, x);
xlabel('time (sec)');
ylabel('output x');

% Now let's do the same thing using Fourier transforms
Y = fft(y);

% w above is angular frequency (radians/sec).  If we
% measure frequencies instead as cycles/sec we need an extra 2pi - thus we
% what we want to solve is
%       X(f) = Y(f) / (alpha + 2 i pi f)

% make frequency axis - remember from Fourier stuff that the frequencies
% are stored in a weird order - from low to high positive frequencies, than
% from high to low negative frequencies.  The two lines below construct
% those frequencies based on a frequency step size of 1/(NumPts * TimeStep)
% - this it the lowest frequency we sample.
Freq(1:NumPts/2) = (0:NumPts/2-1) ./ (NumPts * TimeStep);
Freq(NumPts/2+1:NumPts) = -(NumPts/2 - (0:NumPts/2-1)) ./ (NumPts * TimeStep);

% solve DEQ in frequency domain
X = Y ./ (alpha + i * 2 * pi * Freq);
% put back in time domain
x = ifft(X);
subplot(3, 1, 3);
ylim([0 0.06]);
plot(tme, x);
xlabel('time (sec)');
ylabel('output x');

% these two approaches should give pretty similar answers.  

% The power of this approach is that many differential equations that we
% cannot guess the answer to in the time domain can be solved using the
% Fourier transform.  It is often very useful to have an analytical
% solution because you can identify how the solution depends on particular
% parameters (e.g. alpha above) instead of having to extract that by
% multiple numerical solutions.  



% Numerical solution of differential equations describing 
% transduction cascade in vertebrate rods.  Default
% parameters set for toad rods.

% rate constant for decay of pde activity (units 1/sec)
PDEDecay = 20;
% mean dark pde activity (units 1/sec)
PDEDark = 3;
% shutoff rate for rhodopsin activity (units 1/sec)
RhDecay = 3;
% amplitude of rhodopsin pulse (proportional to flash intensity)
RhAmp = 0.1;
% dark cGMP concentration (units uM)
cGMPDark = 14;
% Ca exchange rate constant (units 1/sec)
ExRate = 0;
% affinity for Ca activation of guanylate cyclase (in uM)
KGC = 0.1;
% cooperativity of action of Ca on cyclase
CoopGC = 2;
% dark Ca concentration (in uM)
CaDark = 0.5;
% constant relating membrane current in pA to cube of cGMP concentration in uM
cGMP2Cur = 8e-3;
% constant relating membrane current to Ca influx
q = ExRate * CaDark / cGMPDark^3;
% maximal rate of cGMP synthesis by cyclase (in uM/sec)
CyclaseMax = (1 + (CaDark / KGC)^CoopGC) * PDEDark * cGMPDark;
% time step for numerical solution of differential equations
TimeStep = 0.001;
% number of points to loop through
NumPnts = 900;
% output arrays for cgmp, pde and ca; start at dark values
clear cgmp;
clear ca;
clear pde;
clear current;
cgmp(1) = cGMPDark;
ca(1) = CaDark;
pde(1) = PDEDark;
current(1) = cGMP2Cur * cGMPDark^3;
% main loop: brute force numerical solution of differential equations
for n = 1:NumPnts
	rh(n) = RhAmp * exp(-n*TimeStep*RhDecay);
	cyclase = CyclaseMax / (1 + (ca(n) / KGC)^CoopGC);
	pde(n+1) = pde(n) + TimeStep * (rh(n) - PDEDecay * (pde(n) - PDEDark));
	cgmp(n+1) = cgmp(n) + TimeStep * (cyclase - pde(n)*cgmp(n));
	ca(n+1) = ca(n) + TimeStep * (q*cgmp(n)^3 - ExRate*ca(n));
	current(n+1) = cGMP2Cur * cgmp(n+1)^3;
end

subplot(2, 1, 1)
plot(current)
subplot(2, 1, 2)
plot(ca)


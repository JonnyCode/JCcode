function [fit] = EPSCfit(coef,x)

% JC 4/19/08

% EPSC(t) = (A(1-e^(-(t-t0)/Taurise))^n)*(Be^(-(t-t0)/Taudecay1)+(1-C)*(e^(-(t-t0)/Taudecay2))
% see Nielsen et al 2004


%fit = coef(1)*((1-exp(-1*(x-coef(2))./coef(3))).^coef(4)).*-(coef(5)*exp(-1*(x-coef(2))./coef(7))+(1-coef(6))*exp(-1*(x-coef(2))./coef(8))) ; %epsc
fit = coef(1)*((1-exp(-1*(x-coef(2))./coef(3))).^coef(4)).*+(coef(5)*exp(-1*(x-coef(2))./coef(7))+(1-coef(6))*exp(-1*(x-coef(2))./coef(8))) ; %ipsc

% fit(1:floor(coef(2))) = 0 ; 

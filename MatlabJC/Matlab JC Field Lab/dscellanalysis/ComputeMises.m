function response = ComputeMises(param, x)

% response = A*exp(kappa*cos(x-u)*pi/180)/exp(kappa) [Tylor]
% A = param(1);
% u = param(2);
% kappa = param(3);

A = param(1);
u = param(2);
kappa = param(3);
response = A.*exp(kappa.*cos(x-u))./exp(kappa);
%response2 = A.*exp(kappa.*cos(x-u)) + m;
end
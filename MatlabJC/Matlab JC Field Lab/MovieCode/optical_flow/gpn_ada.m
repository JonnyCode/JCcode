% Simulate a goal programming network
% Adaptive learning rate strategy has been used
% For details, see
% Van Hulle, M.M.  (1991).  A Goal Programming Network for Linear
% Programming, Bio. Cybern., 65, 243--252.
%
%
% Usage: [V,F,E,U,M,conv_flag] = gpn_ada (D, B, G, ep, dt, U_init, m_init)
%
%		V	Final Network State
%		F	Constraint Satisfaction (F outputs)
%		E	sum(abs(F))
%		D,B	Constraint arrays
%		G	Gain Vector [nc x 1]
%		M	Learning rate trace

function [V,F,E,U,M,conv_flag] = gpn_ada (D, B, G, ep, dt, U_init, m_init)

[nc ns] = size(D);	% # states, # constraints

if (nargin<4)
	ep = 1000;
end
if (nargin<5)
	dt = .001;
end
if (nargin<6)
	U_init = zeros(ns,1);
end
if (nargin<7)
	m_init = 1;
end

conv_flag = 0;

[V F E U M conv_flag] = gpn_ada_mex (D, B, G, ep, dt, U_init, m_init);



function fit = WeberFechner(coef,x)

% weber-fechner function taken from Angueyra and Rieke 2013

% coef(1) = half adapting background
% coef(2) = slope modifier

%fit = 1./(1+(x/coef(1))) ;

fit = 1./(1+(x/coef(1)).^coef(2)) ; % modified weber fechner function

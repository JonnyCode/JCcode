function LJP = liquidJPestimate(val,mob,concPip,concBath, RToverF) 

% this function will estimate the liquid junction potential based on the
% henderson liquid junction potential equation.  taken from Barry and Lynch 1991.
% JC 3/9/12

% val (valence), mob (mobility), and concPip/concBath (concentration in the pipette and bath) should be in vectors
% RToverF (RT/F where R is the gas constant, T is abs temp in kelvins, and
% F is faradays constant) is a constant

% RT/F at 25C is 25.69

Sf = sum(val.*mob.*(concBath-concPip))/sum(val.^2.*mob.*(concBath-concPip)) ;

LJP = RToverF*Sf*log(sum(val.^2.*mob.*concPip)/sum(val.^2.*mob.*concBath)) ;


function Output = NonLinFilter(Input,InputBins,Nonlinearity) ;

% this function will implement the nonlinearity generated in
% NonLinFilterFinder.m on a given input.

% JC 8/14/10


 Output = interp1(InputBins,Nonlinearity,Input,'linear','extrap') ; % assume linearity between bins

    


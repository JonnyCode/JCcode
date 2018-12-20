function SpikePnts = NLPmodel(Stim, ModelParams.Srf, ModelParams.Trf, ModelParams.Nl) ;

% JC 1/14/15

% this function takes a stimulus (Stim) and convolves it with a 2d spatial
% (Srf) and temporal (Trf) receptive field, passes it through a defined
% nonlinearity and generates spike times from a poison process.

% stim should be a 3d (NxNxt) matrix (N spatial and t time points) of
% stimulus intesities at each location and time point.  Srf should be a 2d
% NxN matrix of spatial response sensitivities.  Trf should be a 1d vector
% of temporal sensitivities.  Nl should be a 2 row matrix with the first
% row indicating stimulus intensity and the second row spike rate.


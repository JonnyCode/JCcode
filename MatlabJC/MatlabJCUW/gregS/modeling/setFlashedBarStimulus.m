function M = setFlashedBarStimulus(M,params)
barLocation = params.barLocation;
barFrames = params.barFrames;
nFrames = params.nFrames;
amplitude = params.amplitude;
stimSize = params.stimSize;

S = zeros(stimSize,stimSize,nFrames);
S(:,barLocation,barFrames) = amplitude;

M.stimulus = S;

%Is this the right place for this?
M = applySpatialFilter(M);
M.options = struct;
M.options.preComputedInput = 1;

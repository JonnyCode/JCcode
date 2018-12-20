function Model = setTextureStimulus(Model,params)
%parameters
stimSize = params.stimSize;
textureSigma = params.textureSigma;

%transform
T = TwoD2OneDTransform;
T.rows = stimSize;
T.cols = stimSize;

disp('Adding Stimulus');
s = generateTexture(stimSize,textureSigma);
Model.stimulus = T.invert(s);
disp('done');
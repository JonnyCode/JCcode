% Can we get correlated signals and anticorrelated noise in simple model
% for neural integration?  Consider split field experiment - half of inputs
% see +stm, other half -stm.  Exc and Inh are rectified.  Noise in inputs
% is spatially uncorrelated.  

% FR and JC 
% JC modified 1/11/11

NumTrials = 500;
NumPts = 1000;
Thresh = -.2;          % rectification

% signal - will be + for group1 and - for group2
sig = normrnd(0, 1, NumPts, 1);

clear response_exc response_inh resid_exc resid_inh

for trial = 1:NumTrials
    noise1 = normrnd(0, 1, NumPts, 1);      % group 1 noise
    noise2 = normrnd(0, 1, NumPts, 1);      % group 2 noise

    % create group 1 exc and inh
    exc1 = sig + noise1;
    exc1(find(exc1 < Thresh)) = 0;
    inh1 = -sig - noise1;
    inh1(find(inh1 < Thresh)) = 0;

    % create group 2 exc and inh
    exc2 = -sig + noise2;
    exc2(find(exc2 < Thresh)) = 0;
    inh2 = sig - noise2;
    inh2(find(inh2 < Thresh)) = 0;

    % store integrated responses
    response_exc(trial, :) = exc1 + exc2;  
    response_inh(trial, :) = inh1 + inh2;
end

% how correlated are signals?  
cc_temp=corrcoef(mean(response_exc), mean(response_inh)) ;
cc_signal = cc_temp(1,2) 

for trial = 1:NumTrials
    resid_exc(trial, :) = response_exc(trial, :) - mean(response_exc);
    resid_inh(trial, :) = response_inh(trial, :) - mean(response_inh);
end

% how correlated are residuals?
cc_temp = corrcoef(resid_exc, resid_inh) ;
cc_noise = cc_temp(1,2)
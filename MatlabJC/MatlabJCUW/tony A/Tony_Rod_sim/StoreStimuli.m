function stim = StoreStimuli(EpochCondition,cond)
%function ReturnedEpochCondition = StoreStimuli(EpochCondition,cond)
%
% Creates the stimulus trace from the parameters stored in the
% EpochCondition and then places the vector in
% ReturnedEpochCondition.Stimulus

DEBUG = 0;
Cond = EpochCondition(cond);
dur = FindSearchPara(Cond,'EpochPts');
stim = zeros(1,dur);

stim(:) = Cond.UserInfo.CalibratedMean/1000;

pre = FindSearchPara(Cond,'PrePoints');
stimDur = FindSearchPara(Cond,'StimDur');

stim((pre+1):(pre+stimDur)) = Cond.UserInfo.StimulusAmp;

if (DEBUG)
    view = input('View Stimulus? 0 - no, 1 - yes...');
    if (view)
        h = figure;
        disp(strcat('Displaying: -',EpochCondition(cond).Label{1}));
        plot(1:dur,stim);
        pause;
        close(h);
    end
    clear h;
end
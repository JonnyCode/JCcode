function stm = StoreNoiseStimuli(EpochCondition,file,cond)
% ReturnedEpochCondition = StoreNoiseStimuli(EpochCondition,cond)
%
% Returns a noise stimulus

epoch = EpochCondition(cond).EpochNumbers(1);
[fp, error] = ITCInitializeAnalysis(500000, file);
[stm, error] = ITCReadEpochStm(epoch, 0, fp);

DEBUG = 0;

if (DEBUG)
    view = input('View Stimulus? 0 - no, 1 - yes');
    if (view)
        h = figure;
        disp(strcat('Displaying: ',EpochCondition(cond).Label{1}));
        plot(stm);
        pause;
        close(h);
    end
    clear h;
end
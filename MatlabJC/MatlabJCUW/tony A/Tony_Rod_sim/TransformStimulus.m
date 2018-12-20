function stim = TransformStimulus(EpochCondition,cond)
% function stim = TransformStimuli(EpochCondition,cond)
%
% returns the symetric FFT of the Stimulus vector field of EpochCondition

stim = fft(EpochCondition(cond).StimulusVector...
    (FindSearchPara(EpochCondition(cond),'PrePoints')+1:length(EpochCondition(cond).StimulusVector)));
pstim = fftshift(stim).*conj(fftshift(stim));

DEBUG = 0;
if(DEBUG)
    view = input('View transform? 0 - no, 1 - yes');
    if (view)
        h = figure;    
        disp(strcat('Displaying log of transform: ',EpochCondition(cond).Label{1}));
        plot(FreqAxis(length(pstim)),...
            log(pstim));
        pause;
        try close(h);
        catch
        end
        end
    clear h view;
end
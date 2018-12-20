function xform = TransformData(EpochCondition,cond)
% xform = TransformData(EpochCondition,cond)
% 
% returns the fft for the average of the epochs in EpochCondition(cond)

prepoints = FindSearchPara(EpochCondition(cond),'PrePoints');
xform = fft(EpochCondition(cond).AverageResponse(prepoints+1:length(EpochCondition(cond).AverageResponse)));
pstim = fftshift(xform).*conj(fftshift(xform));

DEBUG = 1;
if(DEBUG)
    view = input('View transform? 0 - no, 1 - yes...');
    if (view)
        h = figure;    
        disp(strcat('Displaying log of transform: ',EpochCondition(cond).Label{1}));
        plot(log(pstim));
        pause;
        
    
        ConvinceMe = input('Convince You? 0 - no, 1 - yes...');
        if(ConvinceMe)
            disp(strcat('Displaying inverse transform: ',EpochCondition(cond).Label{1}));
            plot(ifft(xform));
            pause;
        end
        close(h);
        clear h;
    end
    clear view;
end

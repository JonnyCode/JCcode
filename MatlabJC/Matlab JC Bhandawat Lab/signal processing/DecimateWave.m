function NewWave = DecimateWave(OldWave, StepSize)% Brute force decimation of data in vector.% created FMR ?% edited JC 8/11/08if StepSize~=1 ; % if this is a decimation    NumPts = floor(length(OldWave) / StepSize);    data(1:NumPts) = 0;    for i=1:NumPts        for j=1:StepSize            data(i) = data(i) + OldWave((i-1)*StepSize + j);        end    end    data = data ./ StepSize;    NewWave = data;else    NewWave = OldWave ;     end
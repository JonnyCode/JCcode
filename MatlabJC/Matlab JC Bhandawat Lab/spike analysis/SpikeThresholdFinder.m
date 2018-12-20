function SpikeThreshold = SpikeThresholdFinder(VoltageData, SpikePnts, SI, ThresholdDefinition) ;

% this function will take a matrix of votlageData and cell array of
% spikePnts derived from "SpikeDetection_WC.m" and find the spike threshold of each spike

VoltageData = lowPassFilter(VoltageData,1/SI,2000) ; % lowpass filter voltage data

SpikeTimeSlide = 0.002 ; % seconds past spike point you might fight spike peak
SpikeThreshTimeSlide = .001 ; % seconds before spike peak where spike thresh may be found

SpikePntsSlide = 1/SI * SpikeTimeSlide ; % number of points past spike point you might fight spike peak
SpikeThreshPntsSlide = 1/SI * SpikeThreshTimeSlide ;


if ThresholdDefinition==1 ; % spike Threshold defined as the voltage at the point of the peak of the second derivative preceding a spike

    SpikeThreshold = SpikePnts ; % preallocate for memmory and speed

    for a = 1:length(SpikePnts) % for each spike epoch
        if ~isempty(SpikePnts{a}) ; % if spikes were detected

            VoltageData_diff2 = diff(diff(VoltageData(a,:))) ; % second derivative of voltage data

            for b=1:length(SpikePnts{a}) ; % for every spike
                peak = max(VoltageData(a,SpikePnts{a}(b):SpikePnts{a}(b)+SpikePntsSlide)) ; % spike amplitude at peak
                slidePnts = find(VoltageData(a,SpikePnts{a}(b):SpikePnts{a}(b)+SpikePntsSlide)==peak,1,'first') ; % pnts past detected spike point of peak of spike
                peakPnt = SpikePnts{a}(b)+slidePnts-1 ; % index of peak of spike

                diff_peak = max(VoltageData_diff2(peakPnt-SpikeThreshPntsSlide:peakPnt)) ; % peak of 2nd deriv preceding spike peak
                slidePnts = find(VoltageData_diff2(peakPnt-SpikeThreshPntsSlide:peakPnt)==diff_peak,1,'last') ; % find the local max of the 2nd deriv voltage change that immediatley precedes a spike
                ThreshPnt = peakPnt -SpikeThreshPntsSlide +slidePnts -1 ;

                SpikeThreshold{a}(b) = VoltageData(a,ThreshPnt) ;

            end
        end
    end


elseif ThresholdDefinition==2 ; % spike threshold defined as the highest voltage value acheived that does not lead to a spike or immediatley follow a spike
    
    RefractoryTime = 0.1 ; % seconds after which a spike influences the spike threshold of the voltage that follows
    RefractoryPnts = 1/SI*RefractoryTime ;
    
    maxlocalMax = nan(size(VoltageData)) ; % preallocate for speed
    
    for a = 1:length(SpikePnts) % for each spike epoch
        if ~isempty(SpikePnts{a}) ; % if spikes were detected
            
            for b=1:length(SpikePnts{a}) ; % for every spike     
              
                if b==1 ;
                    [localMaxPnts,localMaxValues] = localMaxFinder(VoltageData(a,1:SpikePnts{a}(b))) ; % find the local maximum before the first spike 
                    maxlocalMax(a,b) = max(localMaxValues) ; % highest local max
                
                elseif b>1 && (SpikePnts{a}(b)-SpikePnts{a}(b-1))>RefractoryPnts ;  
                    [localMaxPnts,localMaxValues] = localMaxFinder(VoltageData(a,SpikePnts{a}(b-1)+RefractoryPnts:SpikePnts{a}(b))) ; % find the local maximum before the first spike     
                    maxlocalMax(a,b) = max(localMaxValues) ; % highest local max
                end
            end
        end
    end
    
    SpikeThreshold = max(maxlocalMax,[],2) ; % highest local max that is not labeled as a spike
        
end

                
classdef FFT1DSpace1DTimeTransform < handle
    %includes a high and low pass cutoff freq (and sampleRate) in time and space
    %y dimension is time
    %input X is 2D, space x time
    properties
        sampleRateTime %samples per second
        sampleRateSpace %points per mm
        cutoffFreqTime_High %in Hz
        cutoffFreqTime_Low
        cutoffFreqSpace_High %in cycles/mm
        cutoffFreqSpace_Low
    end
    methods
        function y = invert(self,x)
            %x_fft = fft2(x);
            [Lspace, Ltime] = size(x);
            
            %space first
            x_fft_space = fft(x,[],1);
            %cutting dimensions we don't need
            FreqStepSize_space = self.sampleRateSpace / Lspace;
            FreqStepSize_time = self.sampleRateTime / Ltime;
            
            %low pass step
            if ~isempty(self.cutoffFreqSpace_High)
                FreqCutoffPts_space = round(self.cutoffFreqSpace_High / FreqStepSize_space);
                x_fft_space(FreqCutoffPts_space:Lspace-FreqCutoffPts_space,:) = 0;
            end
            %high pass step
            if ~isempty(self.cutoffFreqSpace_Low)
                FreqKeepPts_space = round(self.cutoffFreqSpace_Low / FreqStepSize_space);
                x_fft_space(1:FreqKeepPts_space,:) = 0;
                x_fft_space(end-FreqKeepPts_space:end,:) = 0;
            end
            
            %revert
            x = ifft(x_fft_space,[],1);
            %now time fft
            x_fft_time = fft(x,[],2);
            
            %time
            %low pass
            if ~isempty(self.cutoffFreqTime_High)
                FreqCutoffPts_time = round(self.cutoffFreqTime_High / FreqStepSize_time);
                x_fft_time(:,FreqCutoffPts_time:Ltime-FreqCutoffPts_time) = 0;
            end
            
            %high pass
            if ~isempty(self.cutoffFreqTime_Low)
                FreqKeepPts_time = round(self.cutoffFreqTime_Low / FreqStepSize_time);
                x_fft_time(:,1:FreqKeepPts_time) = 0;
                x_fft_time(:,end-FreqKeepPts_time:end) = 0;
            end
            
            %use sparse matrix
            y = sparse(x_fft_time);            
        end
        function y = revert(self,x)
            y = real(ifft(full(x),[],2));
        end
    end
end
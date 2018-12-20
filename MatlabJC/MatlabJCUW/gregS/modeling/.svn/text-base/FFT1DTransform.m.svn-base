classdef FFT1DTransform < handle
    properties
    end
    methods
        function y = invert(self,x)
            y = fft(x);
        end
        function y = revert(self,x)        
             y = real(ifft(x));
        end
    end
end
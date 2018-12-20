classdef TwoD2OneDTransform < handle
    properties
        rows
        cols
    end
    methods
        function y = invert(self,x)
            [self.rows, self.cols] = size(x);
            y = reshape(x,1,self.rows*self.cols);
        end
        function y = revert(self,x)
%             L = self.rows*self.cols;
            %zero padding
%             if length(x)<L
%                 temp = zeros(1,L);
%                 temp(1:length(x)) = x;
%                 x = temp;
%             end
             y = reshape(x,self.rows,self.cols);
        end
    end
end
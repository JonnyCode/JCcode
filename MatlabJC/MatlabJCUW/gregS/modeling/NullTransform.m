classdef NullTransform
    properties
    end
    methods
        function y = invert(self,x)
            y = x;
        end
        function y = revert(self,x)
            y = x;
        end
    end
end
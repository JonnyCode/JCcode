classdef LeafSearchQuery < hgsetget
    properties
        pattern = '';
        fieldnames = {};
        operators = {};
        values = {};
    end
    
    methods        
        function C = LeafSearchQuery()
        end
        
        function C = addCondition(C)
            
            if isempty(C.fieldnames)
                i=1;
                C.pattern = '@1';
            else
                i=length(C.fieldnames) + 1;
                C.pattern = [C.pattern ' && @' num2str(i)];
            end
            
            C.fieldnames{i} = input('Field Name (no quotes): ','s');
            C.operators{i} = input('Operator (no quotes): ','s');
            C.values{i} = input('Value (quotes if a string): ');
            
            displayCondition(C);  
        end
        
        function displayCondition(C)
            
            N = length(C.fieldnames); %number of sub-conditions in search query
            disp('-----------------------');
            for i=1:N
                s = ['@' num2str(i) ': ' C.fieldnames{i} ' ' C.operators{i} ' ' num2str(C.values{i})];
                disp(s);
            end
            
            disp(['Pattern: ' C.pattern]);
            disp('-----------------------');
            
        end
    end
end
    
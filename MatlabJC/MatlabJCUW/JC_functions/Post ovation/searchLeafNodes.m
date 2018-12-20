function selected_leaves = searchLeafNodes(tree,Conditions)

%this function selects leafs as specified by specific parameters. It 
%operates on an object 'Conditions' (created by the class'LeafSearchQuery') 

% GS 1/2/08

N = length(Conditions.fieldnames); %number of sub-conditions in search query
leafMap = tree.leafNodes{1}.splitKeyPath;


%construct query string for each sub-condition
%and place each one in the correct place in Condition.pattern

queryString = Conditions.pattern;
for i=1:N
    if isstr(Conditions.values{i}) %if value is a string
        if strcmp(Conditions.operators{i},'==') %and the condition is equality, do a strcmp
            s = ['strcmp(leafMap.get(' '''' Conditions.fieldnames{i} '''' '),' '''' Conditions.values{i} '''' ')'];
        else
           disp('Error: only equality can be tested for strings'); 
        end
    elseif isnumeric(Conditions.values{i})
        s = ['leafMap.get(' '''' Conditions.fieldnames{i} '''' ')' Conditions.operators{i} num2str(Conditions.values{i})];
    end
    queryString = regexprep(queryString, ['@' num2str(i)], s);
end

disp(queryString)

z = 1; %selected_leaves counter
for i=1:length(tree.leafNodes) %test each leaf
    leafMap = tree.leafNodes{i}.splitKeyPath; 
    if eval(queryString)
        selected_leaves{z} = tree.leafNodes{i};
        z=z+1;
    end
end
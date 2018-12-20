function InterLeafPlot = InterLeafPloter(tree, leafDef) 

%this function will take parameters of splitKeyPath values and plot all
%leaves (mean of epochs within a leaf) that fit that definition

%input: leafDef = a cell with an entry for each splitKey parameter you want
%to search on (leave blank those you don't)

%JC 12/31/08

for a = 1:length(tree.leafNodes) ; % for each leaf in the tree
    h = tree.leafNodes{a}.splitKeyPath ; % this method gets the unique parameters and their values that define leaf 1  

    for b = 1:length(tree.splitKeyPaths) ; % for each parameter   
        if ~isempty(leafDef{b}) ; % if you defined a parameters to search on 
            values = h.get(tree.leafNodes{a}.splitKeyPaths{b}) ; % get the value (this cannot be done using "values(h)" possibly because DC made the map himself?
            
    
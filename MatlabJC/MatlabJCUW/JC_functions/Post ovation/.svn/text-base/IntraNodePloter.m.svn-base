function IntraNodePlot = IntraNodePloter(tree,leafNum) ;

% this function will plot all the means of the epochs that compose a leaf 
% within a single node (i.e. the leaves being compared have the same 
% parent) on a single graph and display the parameters that define each 
% leaf. Additionally, you can click on graph and add key word to epochs of
% the selected leaf 

% input: leafNum = a number of leaf who is a child of the node you want to
% plot

% JC 12/30/08

color = {'b','r','g','k','y','c','m',...
    'b-.','r-.','g-.','k-.','y-.','c-.','m-.',...
    'b:','r:','g:','k:','y:','c:','m:'...
    'b--','r--','g--','k--','y--','c--','m--',...
    'b','r','g','k','y','c','m'...
    'b-.','r-.','g-.','k-.','y-.','c-.','m-.',...
    'b:','r:','g:','k:','y:','c:','m:'...
    'b--','r--','g--','k--','y--','c--','m--'} ;

leaf = tree.leafNodes{leafNum} ;

figure
subplot(4,1,1:3)
hold on
subplot(4,1,4)
hold on

% plot data
subplot(4,1,1:3)
for a = 1:length(leaf.parent.children) ; % for each leaf within the selected node
     LeafMatrix = leaf.parent.children{a}.epochList.responsesByStreamName('Amp_1') ; % get a matrix of epochs within a leaf
     LeafMean{a} = mean(LeafMatrix) ; % get the mean of epochs within the leaf
     NumEpochs(a) = size(LeafMatrix,1) ; % number of epochs per leaf
     samp = leaf.parent.children{a}.epochList.elements{1}.protocolSettings.('acquirino:samplingInterval')/1e6 ; % sample interval in seconds of epoch
     x{a} = [samp:samp:samp*length(LeafMean{a})] ; % x axis
     plot(x{a},LeafMean{a}, color{a})
     xlabel('seconds')
     ylabel('pA or mV')
end
legend(color)

% find parameters
subplot(4,1,4)
PandV = tree.splitKeyPaths ;
for b = 1:length(tree.splitKeyPaths) ; % for each parameter of the map
    for a = 1:length(leaf.parent.children) ; % for each leaf within the selected node
        h = leaf.parent.children{a}.splitKeyPath ; % this method gets the unique parameters and their values that define leaf 1
        v = h.get(tree.splitKeyPaths{b}) ; % values of parameter for leaf a
        PandV{b} = [PandV{b},', ',num2str(v)] ;     
    end
end
PandV{b+1}=['number of epochs:',num2str(NumEpochs)] ; 
text(.1,.5,PandV)  
     
% interactive feature to add key word tag to all elements of a selected leaf
q1 = input('select leaves to add key word tag (y or n)?') ; % check if they want the interactive part
while strcmp(q1,'y') % while they still want to select traces
    display('click leaf trace you want to add key word tag and then hit keyboard')
    pause
    ax = gca ; % get current axis
    cp = get(ax,'CurrentPoint') ; % this gives last point clicked by mouse matrix
    cp = cp(1,1:2) ; % this gives point x,y of click in appropriate units
    for a = 1:length(leaf.parent.children) ; % for each leaf within the selected node  
        ClickDist(a) = min(sqrt(((x{a}-cp(1)).^2+(LeafMean{a}-cp(2)).^2))) ; %distance of click from each plotted leaf
    end
    [minCD,iCD] = min(ClickDist) ; % find the trace closest to the click
    plot(x{iCD},LeafMean{iCD}, 'ro')
    display('selected trace for key word')
    q2 = input('key word tag (in quotes):')
    for a = 1:length(leaf.parent.children{iCD}.epochList.elements) % for each element (epoch) of this leaf
        leaf.parent.children{iCD}.epochList.elements{a}.addKeywordTag(q2) ; % this method adds a key word
    end
q1 = input('select leaves to add key word tag (y or n)?') ; % check if they want to select another trace   
end
   
IntraNodePlot = leafNum ; % becuase a function needs an output (maybe?)





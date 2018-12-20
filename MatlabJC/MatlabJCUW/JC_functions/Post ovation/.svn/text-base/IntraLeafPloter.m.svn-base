function IntraLeafPlot = IntraLeafPloter(tree,leafNum) ;

% this function will plot the epochs that compose a leaf on a single graph,
% note the epoch numbers, and provide the parameters and values that define
% that leaf.

% input: leafNum = leaf number 

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
epochNumStr = cell(1,length(leaf.epochList.elements)) ; % prep matrix for speed
subplot(4,1,1:3)
for a = 1:length(leaf.epochList.elements) ; % for each epoch in the leaf
    data = leaf.epochList.elements{a}.responses.('Amp_1').data ; % data
    samp = leaf.epochList.elements{a}.protocolSettings.('acquirino:samplingInterval')/1e6 ; % sample interval in seconds
    x = [samp:samp:samp*length(data)] ; % x axis
    plot(x,data,color{a}) ;
    xlabel('seconds')
    ylabel('pA or mV')
    epochNum = leaf.epochList.elements{a}.protocolSettings.('acquirino:epochNumber') ; % epoch number
    epochNumStr{a} = num2str(epochNum) ;
end
legend(epochNumStr) ; 

% find parameters
subplot(4,1,4)
h = leaf.splitKeyPath ; % this method gets the unique parameters and their values that define a leaf
params = keys(h) ; % gets the parameters only from the "map", h
for a = 1:length(params) ; % for each parameter of the map
    values = h.get(params{a}) ; % get the value (this cannot be done using "values(params)" possibly because DC made the map himself?
    PandV{a} = [params{a},' : ',num2str(values)] ;    
end
text(.1,.5,PandV)
    
IntraLeafPlot = leafNum ; % becuase a function needs an output (maybe?)








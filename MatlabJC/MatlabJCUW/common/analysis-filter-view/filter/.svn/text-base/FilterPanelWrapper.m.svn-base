function FilterPanelWrapper(tree,figH,doInit)
if doInit %set up for the first time
    UserData = get(figH,'UserData');
    if isfield(UserData,'FilterPanel') && ishandle(UserData.FilterPanel.panelH) %if panel is already there and window is open, just initialize
        %init?
    else
        UserData.FilterPanel = FiltPanel(tree);
        set(figH,'UserData',UserData);
    end
    %wipe out children from previous calls
    ch = get(figH,'children');
    for i=1:length(ch)
        delete(ch(i));
    end
else %run the function
    UserData = get(figH,'UserData');
    UserData.FilterPanel.tree = tree;    
    if ishandle(UserData.FilterPanel.panelH) %if panel found
        %run
        set(figH,'UserData',UserData)
        selection = get(UserData.FilterPanel.handles('selection_toggle'),'Value');
        selection = abs(selection - 1); %invert
        disp('Filter Starting');
        runFilter(UserData.FilterPanel.tree,UserData.FilterPanel.filt.makeQueryString,selection);
        disp('Filter Done');
    else %panel closed, initialize again
        UserData.FilterPanel = FiltPanel(tree);
        set(figH,'UserData',UserData);
    end
end



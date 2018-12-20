function ViewPanelWrapper(tree,figH,doInit)
if doInit %set up for the first time
    UserData = get(figH,'UserData');
    if isfield(UserData,'ViewPanel') && ishandle(UserData.ViewPanel.panelH) %if panel is already there and window is open, just initialize
        %UserData.ViewPanel.init(); %initialize ViewPanel menus (with new tree)
        UserData.ViewPanel.calc.init(); %initialize hidden params and calculation time
    else %make new panel
        UserData.ViewPanel = CalcPanel(tree,figH);
        set(UserData.ViewPanel.panelH,'name','View Panel');
        set(figH,'UserData',UserData);
    end
    %wipe out children from previous calls (ben putting this in gui
    %code)
    ch = get(figH,'children');
    for i=1:length(ch)
        delete(ch(i));
    end
else %run the function
    UserData = get(figH,'UserData');
    UserData.ViewPanel.tree = tree;
    %if ~isfield(UserData.ViewPanel, 'panelH') %need to re-init
    %    UserData.ViewPanel = CalcPanel(tree,figH);
    %    set(UserData.ViewPanel.panelH,'name','View Panel');
    %    set(figH,'UserData',UserData);
    %end
    if ishandle(UserData.ViewPanel.panelH)
        if ~isempty(UserData.ViewPanel.calc) %if there is a calculation loaded
            %this takes a really long time % ok first time
            %UserData.ViewPanel.init(); %initialize ViewPanel menus (with new tree)      
            UserData.ViewPanel.calc.init(); %initialize hidden params and calculation time
            tic;
            tree.accept(UserData.ViewPanel.calc);
            toc;
        end
        
    else %panel was closed, initialize again
        UserData.ViewPanel = CalcPanel(tree,figH);
        set(UserData.ViewPanel.panelH,'name','View Panel');
        set(figH,'UserData',UserData);
    end
end

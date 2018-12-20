function CalcPanelWrapper(tree,figH,doInit)
if doInit %set up for the first time
    UserData = get(figH,'UserData');
    if isfield(UserData,'CalcPanel') && ishandle(UserData.CalcPanel.panelH) %if panel is already there and window is open, just initialize
        %UserData.CalcPanel.init(); %initialize CalcPanel menus (with new tree)
        UserData.CalcPanel.calc.init(); %initialize hidden params and calculation time
    else %make new panel
        UserData.CalcPanel = CalcPanel(tree,figH);
        set(UserData.CalcPanel.panelH,'name','Calculation Panel');
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
    UserData.CalcPanel.tree = tree;
    if ishandle(UserData.CalcPanel.panelH) %panel found
        if ~isempty(UserData.CalcPanel.calc) %if there is a calculation loaded
            %UserData.CalcPanel.init(); %initialize CalcPanel menus (with new tree)
            UserData.CalcPanel.calc.init(); %initialize hidden params and calculation time
            tree.accept(UserData.CalcPanel.calc);
        end
        
    else %panel was closed, initialize again
        UserData.CalcPanel = CalcPanel(tree,figH);
        set(UserData.CalcPanel.panelH,'name','Calculation Panel');
        set(figH,'UserData',UserData);
    end
end

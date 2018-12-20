function results = keywordTool(tree,figH,doInit)
%remember, figH is a panel in epochTreeGUI
results = struct;

%Title Text
pos = [.4 .9 .2 .1];
uicontrol('style','text','tag','keywordTool_text','parent',figH,...
    'string','Keyword Tool','units','normalized','position',pos,...
    'FontSize',18);

%keyword text
pos = [.05 .75 .1 .08];
uicontrol('style','text','tag','keyword_text','parent',figH,...
    'string','Keyword','units','normalized','position',pos,...
    'FontSize',14,'horizontalAlignment','left');

%keyword edit box
pos = [.2 .75 .4 .08];
uicontrol('style','edit','tag','keyword_edit','parent',figH,...
    'string','','units','normalized','position',pos,...
    'FontSize',14,'horizontalAlignment','left');

%add button
pos = [.6 .75 .12 .08];
uicontrol('style','pushbutton','tag','add_button','parent',figH,...
    'string','Add','units','normalized','position',pos,...
    'FontSize',14,'callback',@add_callBack);

%all keywords text
pos = [.05 .65 .15 .08];
uicontrol('style','text','tag','allKeywords_text','parent',figH,...
    'string','All keywords','units','normalized','position',pos,...
    'FontSize',14,'horizontalAlignment','left');

%all keywords popupmenu
pos = [.2 .65 .4 .08];
uicontrol('style','popupmenu','tag','keyword_menu','parent',figH,...
    'string','<select keyword>','units','normalized','position',pos,...
    'FontSize',14,'horizontalAlignment','left');

%remove button
pos = [.6 .65 .12 .08];
uicontrol('style','pushbutton','tag','remove_button','parent',figH,...
    'string','Remove','units','normalized','position',pos,...
    'FontSize',14,'callback',@remove_callBack);

%Keyword field postive for text
pos = [.4 .5 .2 .1];
uicontrol('style','text','tag','addTo_text','parent',figH,...
    'string','Keyword ON for:','units','normalized','position',pos,...
    'FontSize',14);

%button group
pos = [.3 .3 .4 .2];
bgroup = uibuttongroup('tag','bgroup','parent',figH,...
    'title','','units','normalized','position',pos);

%selection radio button
pos = [.05 .5 .4 .3];
uicontrol('style','radiobutton','tag','selection_radioButton','parent',bgroup,...
    'string','selection','units','normalized','position',pos,...
    'FontSize',12);

%selection radio button
pos = [.5 .5 .4 .3];
uicontrol('style','radiobutton','tag','deselection_radioButton','parent',bgroup,...
    'string','deselection','units','normalized','position',pos,...
    'FontSize',12);

%remove all button
pos = [.35 .05 .3 .15];
uicontrol('style','pushbutton','tag','removeAll_button','parent',figH,...
    'string','Remove all keywords','units','normalized','position',pos,...
    'FontSize',14,'callback',@removeAll_callBack);

UserData.handles = makeFigureHandles(figH);
UserData.tree = tree;
set(figH,'UserData',UserData);
if isfield(tree.custom,'keywordMap')
    M = tree.custom.keywordMap;
    set(UserData.handles('keyword_menu'),'String',M.keys);
end

function add_callBack(hObject,eventData)
UserData = get(get(hObject,'parent'),'UserData');
addNodeKeyword(UserData.tree,get(UserData.handles('keyword_edit'),'string'),...
   get(UserData.handles('keyword_edit'),'value'));
M = UserData.tree.custom.keywordMap;
set(UserData.handles('keyword_menu'),'String',M.keys);

function remove_callBack(hObject,eventData)
UserData = get(get(hObject,'parent'),'UserData');
menu_strings = get(UserData.handles('keyword_menu'),'string');
keyword = menu_strings{get(UserData.handles('keyword_menu'),'value')};
removeNodeKeyword(UserData.tree,keyword);
M = UserData.tree.custom.keywordMap;
set(UserData.handles('keyword_menu'),'String',M.keys);
set(UserData.handles('keyword_menu'),'Value',1);

function removeAll_callBack(hObject,eventData)
UserData = get(get(hObject,'parent'),'UserData');
removeAllNodeKeywords(UserData.tree);
set(UserData.handles('keyword_menu'),'String','<select keyword>');
set(UserData.handles('keyword_menu'),'Value',1);


function testGraphicalCheckBox
%Run graphicalCheckBox through some paces.  Don't worry actual graphics.

% should be constructable with no arguments
%%%SU
gBox = graphicalCheckBox();
%%%TS isobject(gBox)
%%%TS isa(gBox, 'graphicalCheckBox')


% should be constructable with axes argument
%%%SU
gBox = graphicalCheckBox(axes);
close all
%%%TS isobject(gBox)
%%%TS isa(gBox, 'graphicalCheckBox')


% should fail to construct on non-axes argument
%%%SU
correctlyFailed = false;
try
    gBox = graphicalCheckBox(5934867);
catch
    correctlyFailed = true;
end
%%%TS correctlyFailed


% invoking boxButtonDownFcn should toggle isChecked property
%   this is the normal click behavior
%%%SU
gBox = graphicalCheckBox(axes);
initialState = gBox.isChecked;
graphicalCheckBox.boxButtonDownFcn(gBox.box, [], gBox);
secondState = gBox.isChecked;
graphicalCheckBox.boxButtonDownFcn(gBox.box, [], gBox);
thirdState = gBox.isChecked;
close all
%%%TS isequal(initialState, thirdState)
%%%TS isequal(initialState, ~secondState)


% boxButtonDownFcn should invoke a callback of the form {@foo}
%%%SU
cbString = 'callback happened';
foo = @(gCheckBox, event) set(gCheckBox.parent, 'UserData', cbString);
gBox = graphicalCheckBox(axes);
gBox.callback = {foo};
graphicalCheckBox.boxButtonDownFcn(gBox.box, [], gBox);
axUserData = get(gBox.parent, 'UserData');
close all
%%%TS strcmp(cbString, axUserData)

% boxButtonDownFcn should invoke a callback of the form {@foo, bar}
%%%SU
foo = @(gCheckBox, event, cbArg) set(gCheckBox.parent, 'UserData', cbArg);
bar = 'callback got argument';
gBox = graphicalCheckBox(axes);
gBox.callback = {foo, bar};
graphicalCheckBox.boxButtonDownFcn(gBox.box, [], gBox);
axUserData = get(gBox.parent, 'UserData');
close all
%%%TS strcmp(bar, axUserData)
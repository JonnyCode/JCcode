% for Andor images
cell = 1;
CellName(cell).name = '071406c2';
cell = cell + 1;
CellName(cell).name = '071406c1';
cell = cell + 1;
CellName(cell).name = '071406c5';

verbose = 0;

for cell = 1:length(CellName)
    PlotConfocalImageG230(CellName(cell).name, verbose);
end


%%
% for confocal images
clear CellName;
cell = 1;
CellName(cell).name = '100808Ac2';
%cell = cell + 1;
%CellName(cell).name = '042308Ac2';
%cell = cell + 1;
%CellName(cell).name = '042308Ac3';

verbose = 1;            % if 1 asks for input as it runs, 0 to run automatically
MaxProjectionFlag = 1;  % should be 1
XYFlag = 1;             % should be 1
LocalFlag = 1;          % local or shared folder for raw and stored images

for cell = 1:length(CellName)
    PlotConfocalImage(CellName(cell).name, MaxProjectionFlag, XYFlag, LocalFlag, verbose);
end
%%


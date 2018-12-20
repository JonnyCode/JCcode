function NeuronStruct2 = NeuronStructGrouper(NeuronStruct,CellTypesArray,paramName1,paramName2)

% This function will allow for two parameter -single group classification of cells in "NeuronStruct" 

% 2016-09-08 JC 

NeuronStruct2 = [] ;

% number of cells
numCells = length(NeuronStruct) ; 

% get CellTypesArray if want all cells
if strcmp(CellTypesArray,'all')
    % get list of unique cell types in NeuronStruct
    CellTypesArray = {} ;
    cnt=1 ;
    for cl = 1:numCells ;
        if ~ismember(NeuronStruct(cl).cellType,CellTypesArray) ; % if this cell type is not in the list
            CellTypesArray{cnt} = NeuronStruct(cl).cellType ;
            cnt = cnt+1 ;
        end
    end
end

% number of cell types 
numCellTypes = length(CellTypesArray) ;

% set color of cell points
ColorMatOrig = [0, 0.4470, 0.7410
    0.8500, 0.3250, 0.0980
    0.9290, 0.6940, 0.1250
    0.4940, 0.1840, 0.5560
    0.4660, 0.6740, 0.1880
    0.3010, 0.7450, 0.9330
    0.6350, 0.0780, 0.1840] ;

ColorMat = ColorMatOrig ;
for a=2:ceil(numCellTypes/size(ColorMat,1)) ;
    ColorMat = [ColorMat; ColorMatOrig/a] % make them darker to expand color space 
end

% get param values
for cl=1:numCells ; % for each cell
    v1(cl) = eval(['NeuronStruct(',num2str(cl),').',paramName1]) ; % value in paramName1
    v2(cl) = eval(['NeuronStruct(',num2str(cl),').',paramName2]) ; % value in paramName2
end
    
% plot histograms and points
figure
subplot(2,2,1)
[h, hx] = hist(v2,100) ;
plot(h,hx)
xlabel('number cells')
ylabel(paramName2)

subplot(2,2,4)
[h, hx] = hist(v1,100) ;
plot(hx,h)
xlabel(paramName1)
ylabel('number cells')

subplot(2,2,2)% plot the points
for ctype = 1:length(CellTypesArray) ; % for each cell type 
    for cl=1:numCells ; % for each cell
        if strcmp(NeuronStruct(cl).cellType,CellTypesArray{ctype}) ; % if its this cell type
            plot(v1(cl),v2(cl),'color',ColorMat(ctype,:),'Marker','*')
            hold on
        end
    end
end
xlabel(paramName1)
ylabel(paramName2)

subplot(2,2,3) % list the cell types names and colors
for ctype = 1:length(CellTypesArray) ; % for each cell type  
    text(0.1,ctype/length(CellTypesArray),CellTypesArray{ctype},'color',ColorMat(ctype,:),'units','normalize')
end
    
% select and rename
subplot(2,2,2)
[x, y] = ginput;
plot(x, y);
IN = inpolygon(v1, v2, x, y);
I = find(IN == 1);
plot(v1(I), v2(I), 'ro')

InputTxt = input('rename cell type (or press enter to not rename):','s') ;
if ~isempty(InputTxt)
    NeuronStruct2 = NeuronStruct ;
    for a=1:length(I)
        NeuronStruct2(I(a)).cellType = InputTxt ;
    end
end


function layer = NullFilter(layer,params)
units = layer.subunits;

disp(['Adding null filter to ' num2str(length(units)) ' subunits...']);
for i=1:length(units);        
    units(i).filter = 1;
end
disp('done');

layer.subunits = units;


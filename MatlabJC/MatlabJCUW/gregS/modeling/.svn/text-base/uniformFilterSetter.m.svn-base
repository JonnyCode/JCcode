function layer = uniformFilterSetter(layer,params)
filter = params.filter; %function handle
units = layer.subunits;

disp(['Adding filter to ' num2str(length(units)) ' subunits...']);
for i=1:length(units);        
    units(i).filter = filter;
end
disp('done');

layer.subunits = units;

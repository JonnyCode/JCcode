function layer = uniformNLSetter(layer,params)
outputNL = params.outputNL; %function handle
units = layer.subunits;

disp(['Adding nonlinearity to ' num2str(length(units)) ' subunits...']);
for i=1:length(units);        
    units(i).outputNL = outputNL;
end
disp('done');

layer.subunits = units;

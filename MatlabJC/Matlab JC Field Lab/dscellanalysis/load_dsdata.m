function[datarun] = load_dsdata(data_path, path2, dg, dg_path, sta)

%Load datarun structure to be analyzed       %sravi 12-17-2012

%Inputs: 
            %data_path: Directory that has all the data folders ('/Analysis/sravi/2012-10-15-0/data000-3600-7200s/datamaps-sr-model/')
            
            %path2: Data directory that needs to be analyzed ('data001-map/data001-map')
            
            %dg: 1 if loading drifting grating data, 0 otherwise
            
            %dg_path: path in which drifting grating stimuli is stored ('/stimuli/s02')
            
            %sta: 1 if datarun has stas and rfs that need to be loaded, 0 otherwise
            
%Output: Datarun structure
 
% Identify what to load
opt = struct('verbose',1,'load_params',1,'load_neurons',1,'load_ei',1);
 %opt = struct('verbose',1,'load_neurons',1);
datarun = load_data([data_path,path2], opt);
 
% parse drifting grating stimulus into matlab

if(dg == 1)
    datarun.names.stimulus_path = [data_path,dg_path];
    datarun = load_stim(datarun);
    %datarun = load_stim(datarun, 'user_defined_trigger_interval', 10);
end
 
%Get STA and receptive field

if(sta == 1)
    datarun = load_sta(datarun, 'load_sta', 'all');
    marks_params.thresh = 4.0;
    datarun = get_sta_summaries(datarun, 'all', 'marks_params', marks_params);
    datarun = get_sta_fits_from_vision(datarun, 'all');
end

end
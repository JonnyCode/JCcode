function datarun = convert_stim_MatlabDG(datarun)
stim.params.SPATIAL_PERIOD = datarun.stimulus.params.spatial_period;
stim.params.TEMPORAL_PERIOD = datarun.stimulus.params.temporal_period;
stim.params.DIRECTION = datarun.stimulus.params.direction;

for i = 1:length(datarun.stimulus.trials)
    stim.trials(i).SPATIAL_PERIOD = datarun.stimulus.trials(i).spatial_period;
    stim.trials(i).TEMPORAL_PERIOD = datarun.stimulus.trials(i).temporal_period;
    stim.trials(i).DIRECTION = datarun.stimulus.trials(i).direction;
end

for i = 1:length(datarun.stimulus.combinations)
    stim.combinations(i).SPATIAL_PERIOD = datarun.stimulus.combinations(i).spatial_period;
    stim.combinations(i).TEMPORAL_PERIOD = datarun.stimulus.combinations(i).temporal_period;
    stim.combinations(i).DIRECTION = datarun.stimulus.combinations(i).direction;
end
datarun.stimulus.params = stim.params;
datarun.stimulus.trials = stim.trials;
datarun.stimulus.combinations = stim.combinations;
end

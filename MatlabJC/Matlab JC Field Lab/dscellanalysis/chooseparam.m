function[mag dsindex magmax magave angle rho theta num NumSpikesCell, StimComb] = chooseparam(datarun)

[NumSpikesCell, StimComb] = get_spikescellstim(datarun, datarun.cell_ids, 0);
%[NumSpikesCell, StimComb] = get_avemaxfiringrate(datarun);
%[NumSpikesCell, StimComb] = get_maxavefiringrate(datarun);
%[NumSpikesCell, StimComb] = get_meanISI(datarun);
%[NumSpikesCell, StimComb] = get_medianISI(datarun);

[mag dsindex magmax magave angle rho theta num] = dscellanalysis(NumSpikesCell, StimComb);

end
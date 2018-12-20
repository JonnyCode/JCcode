
clear all
close all

% seizure anaylsis script to run seizure analysis function
[mean_isi_data,spk_rastor4_zero_data]=seizure_analysis('i070612cData2','ro');% run seizure analysis

% nearest neihbor analysis function call
[number_exacts, number_matches, percent_match_ms, percent_match_quartms]=nearest_neighbor(spk_rastor4_zero_data1,spk_rastor4_zero_data2) ;

% create cross correlation of recorded pair for each sweep
cc = xcov(spk_rastor4_zero_data1,spk_rastor4_zero_data2,'coef') ;

lags = 1:length(cc) ;                  % creates x-axis for lags
lags = lags - (length(cc)+1)/2 ;
lags =  lags*(1/20) ;                  % make lags in ms (1/samplerate*1000=1/20)

figure(4)
plot(lags,cc)

midpoints=find(lags>-5 & lags<5); %find lags between -5 and 5 ms

midsec = cc(midpoints);        % define section to fit with sinwave
midsec = midsec - mean(midsec);    % subtract mean

midlags = lags(midpoints);      % define section of lag to be fit with sinwave

figure(5),plot(midlags,midsec)

cftool  % use to fit with cosine function (fit= A*cos(B*midlags+C) , C = phase

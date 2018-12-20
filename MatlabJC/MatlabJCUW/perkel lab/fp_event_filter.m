% fp_event_filter
% Jon Cafaro
% 12/14/06

%      This script will take detected field potential events from a Clampfit results
% sheet(cut/paste into an Excel file), make sure each event is only detected one 
% time and extract only those events above a determined absolute amplitude
% and instantaneous frequency.  Additionally it will calculate new interterevent
% intervals and instantaneous frequencies.  It will also plot the following
% graphs: 
%       fig (1) trace number v. instantaneous frequency, v. peak amplitude, v.
%       rise slope (10%-90%), v. decay slope (10%-90%).(with averages included)
%       fig (2) trace number v. event time of peak
%       fig (3) same as figure 1 but ploting absolute value
%       fig (4) ratio of catagory 1: catagory 2 events
%       fig (5) number of events in each sweep

% columns: trace number=1, peak amplitude=8, time of peak=10, instantaneous frequency=35, interevent interval=34, rise
% slope=24, decay slope=26


clear, close all
%opens excel file and sets it to be new variable 'raw_events'
cd C:\MATLAB7\rawevents                                                    % change directory to where excel files are stored
file=input('enter the excel spreadsheet file name with the event data in quotes:') %user prompt
raw_events=xlsread(file);                                                  % reads .xls file from user prompt and makes raw_events


% Allow user to enter the "peak amplitude" threshold to create a variable
%"amp_thresh"
amp_thresh=input('filter out all events with an absolute peak amplitude below:')

%   Allow user to enter the "instantanous frequency" threshold to create a variable "inst_freq_thresh" 
inst_freq_thresh=input('filter out all events with instantaneous frequency above (Hz):')


%   Identify events with the same "sweep number" and "time of peak" and remove whichever event has the
%smaller absolute "peak amplitude" creating a new variable "raw_events2".
raw_events2=raw_events;                                                    % sets new variable
j=2;                                                                       % starts j at 2
while j<=size(raw_events2,1);                                              % makes script stop when j has reached end of rows
   if raw_events2(j,1:9:10)==raw_events2(j-1,1:9:10);                      % if column1 and 10 are equal then...   
       if abs(raw_events2(j,8))>abs(raw_events2(j-1,8));                   %   check to see which has larger peak amp
           raw_events2(j-1,:)=[];                                          %   delete row with smaller peak amp   
       else raw_events2(j,:)=[];                                           %   deleter row with smaller peak amp
       j=j;                                                                %   if a row is deleted it moves to next row without skipping 
       end;
   else j=j+1;                                                             %   if a row is not deleted it moves to next row
   end;
end;


% Remove all events with "peak amplitude" less than "amp_thresh" creating a new variable "amp_thresh_events".
amp_thresh_events=raw_events2;                                             % sets new variable
k=1;                                                                       % starts k at 1
while k<=size(amp_thresh_events,1);                                        % keeps k<=number of rows in matrix
    if abs(amp_thresh_events(k,8))<amp_thresh;                             % if the peak amp< amp threshold previously set then...
        amp_thresh_events(k,:)=[];                                                 %delelte that row event
        k=k;                                                                       % moves to next row if row is deleted without skipping
    else k=k+1;                                                            % moves to next row if no row is deleted
    end;
end;


% Calculate new "interevent interval" and "instantaneous frequencies" from
% new events creating a new variable "ampthresh_events2".
amp_thresh_events2=amp_thresh_events;

%calculate new interevent interval and put in new column (34)
amp_thresh_events2(1,34)=NaN;                                              % makes the first event iei NaN
for p=2:size(amp_thresh_events2,1);                                        % sets p=2 to number of rows
    if amp_thresh_events2(p,1)==amp_thresh_events2(p-1,1);                 % if the trace number is the same as the previous event trace number then...
        amp_thresh_events2(p,34)= amp_thresh_events2(p,10)-amp_thresh_events2(p-1,10);  %calculate new interevent interval (iei)
    else amp_thresh_events2(p,34)=NaN;                                     % makes sure iei of first event in trace = NaN
    end;
end;

%calculate new instantaneous frequency and put new column (35)
amp_thresh_events2(1,35)=NaN;                                              % makes the first event if =NaN
for p=2:size(amp_thresh_events2,1);                                        % sets p=2-number of rows
    if amp_thresh_events2(p,1)==amp_thresh_events2(p-1,1);                 % if the trace number is the same as the previous event then...
        amp_thresh_events2(p,35)= (1./(amp_thresh_events2(p,10)-amp_thresh_events2(p-1,10)))*1000.;   %calculates new inst frequ between event and previous event
    else amp_thresh_events2(p,35)=NaN;                                     % makes sure no inst freq will be calc from first event 
    end;
    
end;


%   Identify any events with "instanteous frequency" greater than
% "inst_freq_thresh".  Compare these events to the events directly
% previous and remove the event with the smaller absolute "peak
% amplitude" creating a new variable "freqthresh_events"
freq_thresh_events=amp_thresh_events2;                                      % creates new variable
q=2;                                                                        % sets q=2
while q<=size(freq_thresh_events,1);                                        % keeps q<= row number
    if freq_thresh_events(q,35)>inst_freq_thresh;                           % if inst freq> threshold then
        if abs(freq_thresh_events(q,8))>abs(freq_thresh_events(q-1,8));     %    id the event (q or q-1)with the smaller abs(peak amp) 
            freq_thresh_events(q-1,:)=[];                                   %    delte event with smaller peak amp
        else freq_thresh_events(q,:)=[];                                    %    delete event with smaller peak amp
        q=q;                                                                % move to next row if row is deleted without skipping
        end;
    else q=q+1;                                                             % move to next row if row is not delteted
    end;
end;
        

%Calculate new "interevent interval" and "instantaneous frequencies from
%new events creating a new varaible "filtered_events".
 filtered_events=freq_thresh_events;

%calculate new interevent interval and replace in column (34)
filtered_events(1,34)=NaN;                                                 % sets the first event iei=NaN
filtered_events(1,35)=NaN;                                                 % sets the first event if=NaN
for r=2:size(filtered_events,1);                                           % sets r=1-the number of rows
    if filtered_events(r,1)==filtered_events(r-1,1);                       % makes sure iei is only calc within trace
        filtered_events(r,34)= filtered_events(r,10)-filtered_events(r-1,10); % calculates interval between event and previous event
    else filtered_events(r,34)=NaN;                                           % makes sure iei=0 for the first event of each trace
    end;

%calculate new instantaneous frequency and put new column (35)
    if filtered_events(r,1)==filtered_events(r-1,1);                       %makes sure if are only calculated within a trace
        filtered_events(r,35)= (1./(filtered_events(r,10)-filtered_events(r-1,10)))*1000.;   %calculates ins freq from 1/iei x 1000
    else filtered_events(r,35)=NaN;                                          %makes sure inst freq for first event in trace=NaN
    end;
end;


% calculate averages for each swp number (column 1) of each column a
% creates new variable 'av_filtered_events'
for s=1:max(filtered_events(:,1));                                         % sets s= 1-the number of traces
    trc{s}=find(filtered_events(:,1)==s);                                  % creates trc(x) = all data rows from trace x
    av_filtered_events(s,:)=nanmean(filtered_events(trc{s},2:end),1);      % takes mean of all columns only within a trace x and disregards NaNs
end;
t=[1:max(filtered_events(:,1))]';                                          % creates vector t = sweep numbers
av_filtered_events=[t,av_filtered_events];                                 % vector t to begining of average matrix

% creates columns of event times for each sweep and puts into variable
% 'times_filtered_events'
clear trc;
times_filtered_events=NaN(50,50);                                          % creates a matrix of NaN to put data into
for s=1:max(filtered_events(:,1));                                         % sets s= 1-the number of traces
    trc{s}=find(filtered_events(:,1)==s);                                  % creates trc(x) = all data columns from trace x
    sz=size(filtered_events(trc{s},10),1);                                 % creates a value for the number of rows in the column
    times_filtered_events(1:sz,s)=filtered_events(trc{s},10);              % puts the events times into the NaN matrix
end;                                                           


%calculates duration of events (from peak time of first to last event)
clear trc
for s=1:max(filtered_events(:,1));                                         % see text above
    trc{s}=find(filtered_events(:,1)==s);                                  % see text above
    duration_filtered_events(s,1)=filtered_events(max(trc{s}),10)-filtered_events(min(trc{s}),10);       %subtracts the time of peak of the first event from time of peak of the last event
end
duration_filtered_events=[t,duration_filtered_events];                     % see text above
zip=find(duration_filtered_events(:,2)==0);                                % finds any with zero for duration
duration_filtered_events(zip,2)=NaN ;                                      % turns them into NaN
    

%calculates average duration
av_duration_filtered_events= nanmean(duration_filtered_events);


% replaces all values with absolute values and creates new vairable
% 'absolute_filtered_events'
absolute_filtered_events=abs(filtered_events);

% calculate averages for each trace number in abolute_filtered_events (column 1) 
% of each column and creates new variable 'av_absoulte_filtered_events'
clear trc
for u=1:max(absolute_filtered_events(:,1));                                % sets u= 1-the number of traces
    trc{u}=find(absolute_filtered_events(:,1)==u);                         % creates trc(x) = all data columns from sweep x
    av_absolute_filtered_events(u,:)=nanmean(absolute_filtered_events(trc{u},2:end),1);  % takes mean of all columns only within a trace x and disregards NaNs
end;
v=[1:max(absolute_filtered_events(:,1))]';                                 % creates vector v = trace numbers
av_absolute_filtered_events=[v,av_absolute_filtered_events];               % combines v and av_ab...


%    this finds all events in a trace that are type 1 and type 2 and creates a
% ratio number of type 1/type 2 in the second column of 'cat1v2'.  the first
% column is the trace number.  Also, it tells how many total events per
% sweep.
for w=1:max(filtered_events(:,1));                                         % sets w= 1-the number of traces
    cat1=find((filtered_events(:,1)==w)&(filtered_events(:,3)==1));        % finds index where trace(w) is type 1 event
    cat2=find((filtered_events(:,1)==w)&(filtered_events(:,3)==2));        % finds index where trace(w) is type 2 event
    x=length(cat1);                                                        % calulates the number of type1 events in a trace
    y=length(cat2);                                                        % calcs the number of type 2 events in a trace
    num_filtered_events(w,1)=x+y;                                          % calcs the number of total events per trace
    cat1v2(w,1)=x/y;                                                       % creates a matrix with number type1/number type2 
end;
cat1v2=[v,cat1v2];                                                          % adds a vector for trace number
num_filtered_events=[v,num_filtered_events];                                % adds a vector for trace number
% NOTE- this script can create inf values when there are no type 2 event
% and MATLAB gies a 'dividing by zero' warning.


% gets mean and stdev of av_filtered_events - and creates av2_filtered_events, 
% and stdev2_filtered_events.  It is important to remember this is
% not the same as getting average of all events, instead it is calculating
% the average of the average stats of each sweep.
av2_filtered_events(1,:)=nanmean(av_filtered_events,1);                    % calculates means, puts in row 1
stdev2_filtered_events(1,:)=std(av_filtered_events,1);                     % calcs stand deviation puts in row 1


% gets mean and stdev of av_absolute_filtered_events- remmebber as above
av2_absolute_filtered_events(1,:)=nanmean(av_absolute_filtered_events,1);  % calculates means puts in row 1
stdev2_absolute_filtered_events(1,:)=std(av_absolute_filtered_events,1);   % calcs stand deviation puts in row 2

% calculates the change in number events from sweep #1 vs sweep #4 
% (number of events in sweep#1-number in sweep#4)/number sweep#1
if size(num_filtered_events,1)<4;                                          % if there is no events in the fourth sweep than
    Chng_num_filtered_events=1;                                            % make the percent decay=1
else Chng_num_filtered_events=((num_filtered_events(1,2)-num_filtered_events(4,2))/num_filtered_events(1,2));  % otherwise calculate as mentioned above
end;

% average number of events per sweep
av_num_filtered_events=mean(num_filtered_events);


%calculates latency of the first event in a sweep 
clear trc
for s=1:max(filtered_events(:,1));                                         % 
    trc{s}=find(filtered_events(:,1)==s);                                  % 
    latency_filtered_events(s,1)=filtered_events(min(trc{s}),10);          %
end
latency_filtered_events=[t,latency_filtered_events];                       % 

%calculates average duration
av_latency_filtered_events= nanmean(latency_filtered_events);




%PLOTTING DATA
%plotting from filteres_events
% create a figure with four graphs for trace v. individual peak amp,
% instfreq, rise and decay slopes
figure(1)                                                                  % creates fig1
subplot(2,2,1)                                                             % creates space for sets of axis on fig1

%plots trace number v. peak amplitude
plot(filtered_events(:,1),filtered_events(:,8),'bd')                       % plots
xlabel('trace number')                                                     % x label
ylabel('peak amplitude')                                                   % y label
hold on                                                                    % allows mean data to be graphed on top

%plots trace number v. instantaneous frequency
subplot(2,2,2)                                                             % puts next graph in position 2
plot(filtered_events(:,1),filtered_events(:,35),'bo')                      % plots
xlabel('trace number')                                                     % xlabel
ylabel('instantaneous frequency')                                          % y label
hold on                                                                    % allows mean data to be graphed on top

%plots trace number v. decay slope
subplot(2,2,3)                                                             % puts next graph  in postion 3
plot(filtered_events(:,1),filtered_events(:,26),'b*')                      % plots
xlabel('trace number')                                                     % x label
ylabel('decay slope')                                                      % y label
hold on                                                                    % allows mean data to be graphed on top

%plots trace number v. rise slope
subplot(2,2,4)                                                             % puts next graph in position 4
plot(filtered_events(:,1),filtered_events(:,24),'b+')                      % plots
xlabel('trace number')                                                     % x label
ylabel('rise slope')                                                       % y label
hold on                                                                    % allows mean data to be graphed on top

%plotting from av_filtered_events
%plots trace number v. averages of peak amp, inst freq, rise and decay slp
subplot(2,2,1)                                                             % creates 2x2 graph, puts graph in position 1

%plots trace number v. mean peak amplitude
plot(av_filtered_events(:,1),av_filtered_events(:,8),'rd')                 % plots
xlabel('trace number')                                                     % x label
ylabel('peak amplitude')                                                   % y label
title(file)                                                                % makes title excel file name

%plots trace number v. mean instantaneous frequency
subplot(2,2,2)                                                             % puts next graph in position 2
plot(av_filtered_events(:,1),av_filtered_events(:,35),'ro')                % plots
xlabel('trace number')                                                     % x label
ylabel('instantaneous frequency')                                          % y label
title(file)                                                                % makes title excel file name

%plots trace number v. mean decay slope
subplot(2,2,3)                                                             % puts next graph in position 3
plot(av_filtered_events(:,1),av_filtered_events(:,26),'r*')                % plots
xlabel('trace number')                                                     % x label
ylabel('decay slope')                                                      % y label
title(file)                                                                % makes title excel file name

%plots trace number v. mean rise slope
subplot(2,2,4)                                                             % puts next graph in position 4
plot(av_filtered_events(:,1),av_filtered_events(:,24),'r+')                % plots
xlabel('trace number')                                                     % xlabel
ylabel('rise slope')                                                       % y label
title(file)                                                                % makes title excel file name


%plotting from filtered_events
%plots trace number v. time of peak
figure(2)                                                                  % creates fig2
plot(filtered_events(:,10),filtered_events(:,1),'.')                       % plots
xlabel('time of event peak')                                               % xlabel
ylabel('trace number')                                                     % ylabel
title(file)                                                                % makes title excel filename



%plotting from absolute_filteres_events
% create a figure with four graphs for trace v. individ event peak amp,
% instfreq, rise and decay slope
figure(3)                                                                  % creates fig1
subplot(2,2,1)                                                             % creates space for sets of axis on fig1

%plots trace number v. peak amplitude
plot(absolute_filtered_events(:,1),absolute_filtered_events(:,8),'gd')     % plots
xlabel('trace number')                                                     % x label
ylabel('peak amplitude')                                                   % y label
hold on                                                                    % allows mean data to be graphed on top

%plots trace number v. instantaneous frequency
subplot(2,2,2)                                                             % puts next graph in position 2
plot(absolute_filtered_events(:,1),absolute_filtered_events(:,35),'go')    % plots
xlabel('trace number')                                                     % xlabel
ylabel('instantaneous frequency')                                          % y label
hold on                                                                    % allows mean data to be graphed on top

%plots trace number v. decay slope
subplot(2,2,3)                                                             % puts next graph  in postion 3
plot(absolute_filtered_events(:,1),absolute_filtered_events(:,26),'g*')    % plots
xlabel('trace number')                                                     % x label
ylabel('decay slope')                                                      % y label
hold on                                                                    % allows mean data to be graphed on top

%plots trace number v. rise slope
subplot(2,2,4)                                                             % puts next graph in position 4
plot(absolute_filtered_events(:,1),absolute_filtered_events(:,24),'g+')    % plot
xlabel('trace number')                                                     % x label
ylabel('rise slope')                                                       % y label
hold on                                                                    % allows mean data to be graphed on top

%plotting from av_absolute_filtered_events
%plots trace number v. averages of peak amp, inst freq, rise and decay slp
subplot(2,2,1)                                                             % creates 2x2 graph, puts graph in position 1

%plots trace number v. mean peak amplitude
plot(av_absolute_filtered_events(:,1),av_absolute_filtered_events(:,8),'rd') % plots
xlabel('trace number')                                                     % x label
ylabel('peak amplitude')                                                   % y label
title(file)                                                                % makes title excel file name

%plots trace number v. mean instantaneous frequency
subplot(2,2,2)                                                             % puts next graph in position 2
plot(av_absolute_filtered_events(:,1),av_absolute_filtered_events(:,35),'ro') % plots
xlabel('trace number')                                                     % x label
ylabel('instantaneous frequency')                                          % y label
title(file)                                                                % makes title excel file name

%plots trace number v. mean decay slope
subplot(2,2,3)                                                             % puts next graph in position 3
plot(av_absolute_filtered_events(:,1),av_absolute_filtered_events(:,26),'r*') % plots
xlabel('trace number')                                                     % x label
ylabel('decay slope')                                                      % y label
title(file)                                                                % makes title excel file name

%plots trace number v. mean rise slope
subplot(2,2,4)                                                             % puts next graph in position 4
plot(av_absolute_filtered_events(:,1),av_absolute_filtered_events(:,24),'r+') % plots
xlabel('trace number')                                                     % xlabel
ylabel('rise slope')                                                       % y label
title(file)                                                                % makes title excel file name

set(gcf,'Name','absolute values')                                          % sets name of figure 3


%plots from cat1v2
%plots trace number versus type1 v type 2 ratio
figure(4)                                                                  % creates figure 4
plot(cat1v2(:,1),cat1v2(:,2),'ro-')                                        % plots
xlabel('trace number')                                                     % x label
ylabel('# cat1/# cat2')                                                    % y label
title(file)                                                                % makes title excel filename

%plots from num_filtered_events
%plots the number of events per trace
figure(5)                                                                  % creates fig 5
plot(num_filtered_events(:,1),num_filtered_events(:,2),'ko-')
xlabel('trace number')                                                     %xlabel
ylabel('number of events')                                                 % y label
title(file)                                                                % makes title excel filename


%plots from filtered_events
%plots a histogram of event times
figure (6)                                                                 % creates fig 6
binner=[0:50:50000]                                                        % creates a vector 0 to 50000 skipping by 50
hist(filtered_events(:,10),binner)                                         % creates a histogram of event times binned every 50ms 
xlabel('event time')                                                       % x label
ylabel ('number of events')                                                % y label

%plots from latency_filtered_events
% plots latency v. trace number
figure(7)
plot(latency_filtered_events(:,1),latency_filtered_events(:,2),'go')
xlabel ('sweep number')
ylabel ('latency of first event')

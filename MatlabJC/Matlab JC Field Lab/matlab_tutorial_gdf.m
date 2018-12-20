% The purpose of this tutorial is to introduce you to how to load data and 
% perform some basic kinds of analyses.

% To begin with we will create a datarun structure. This is a structure
% that will be used in most of the analyses of multielectrode array (MEA)
% data in the lab.

% First you should look at the documentation for load_data.  Recall that
% you must AFP mount the server (brahms) in order for matlab to find these
% functions.

help load_data

% The simplest way to use load_data (there are many ways to call this
% function) is as follows:

datarun = load_data('/Volumes/Berlioz/Analysis/2014-09-23-0/data000/data000');
%%%%%%%%datarun = load_data('/Volumes/lab/Experiments/Array/Analysis/2013-11-27-0/data000/data000');

% This function called in this way, does not actually load any data, it
% simply creates a data run MATLAB structure that contains the information to load
% specific pieces of analyzed data.  For more information on matlab
% structures see:
% http://www.mathworks.com/help/matlab/structures.html

% You should note, that I am using data
% somewhat loosely here.  The data you will be loading is pre-analyzed,
% meaning that you will not be loading "raw data" -- voltage as a function
% of time.  Instead you are loading spike times, stas, electrical images
% (EIs), etc.  The data run structure is meant to be an entity in matlab
% that can contain these different kinds of analyzed data so that you can
% further analyze the data. For example, let's say that you want to know
% that average receptive field diameter of all of the ON cells in a given
% recording.  While the Vision software package will calculate the STAs for
% you, it will not calculate this particular quantity.  So you can would
% need to bring the STAs into matlab, gather up all the ON cells, figure
% out the receptive field diameter of each cell and then calculate the
% average.  OK, getting back to our task of learning how to load data into
% matlab. 

% The structure "datarun" contains one element -- another structure called
% names

datarun.names


% To load the spike times, you need will use the function load_neurons,
% because the neurons file from Vision contains this infromation:

help load_neurons

datarun = load_neurons(datarun);

% Now datarun contains several other items besides "names".  It contains,
% "spikes" and "triggers" and other things too.

datarun

% Let's consider further datarun.spikes. 

datarun.spikes

% datarun.spikes is a cell array:
% (http://www.mathworks.com/help/matlab/matlab_prog/what-is-a-cell-array.html)
% It contains all the spikes for all the neurons in this recording.
% If you run...

datarun.spikes{1}

% ...you will get a long list of numbers.  This is the number of spikes
% associate with the first ganglion cell in the datarun structure. To see
% how many spikes there are for this neuron you can type

number_of_spikes = length(datarun.spikes{1})

% To get the number of spikes for the 100th neuron, what do you need to
% type?

% For the first exercise, you should do the following.  For the first
% neuron, plot a histogram of the spikes in 1 second bins.  Plot the 
% histogram using 5s bins, 0.1 s bins, 0.5 s bins, and 10 s bins.  How does
% the histogram change with the bin size?

%%
% Let's now return to the datarun. Look at datarun.triggers.

datarun.triggers

% This is a long list of numbers (a vector).  These are the times when the
% computer that generates the visual stimulus sent a TTL pulse (a trigger)
% to the data acquisition (DAQ) computer.  Depending on the stimulus that was
% presented, the meaning of these triggers can be somewhat different, but
% in general, their purpose to align the timing of visual stimuli with the
% spike times.

% At this point, you should look at the experimental notebook for this data
% set.  You can find it at /Volumes/lab/Experiments/Notebooks/2013-11-27-0/

% This is a recording from a dopamine receptor type 1 knock out (KO)
% mouse.  In data000, the first bit of recorded data from this animal, the
% stimulus was pretty simple: 50% contrast flashes on top of a mean
% background.  These flashes were full-field.  You should think about what
% 50% contrast and full-field mean in this context.

% This is the time of the first flash.
datarun.triggers(1)

% There is a trigger for when the flash goes on AND when it goes off, so
% the flash goes off at

datarun.triggers(2)

% The flash presented by a video display that has a frame rate of 60 Hz.
% How many frames of the display were used to present the flash?

% Notice that the notebook entry says that the stimulus started at NDF 2.
% This means that the mean light intensity was attenuated by a factor of
% 100 (2 log units).  At ~50s we switched the NDF to 0, meaning the mean
% light intensity increased by 100 fold.  The purpose of this experiment
% was to see how much the response to the 50% contrast flash would change
% after switching from rod mediated responses (NDF 2) to cone mediated
% responses (NDF 0).  

% Your exercise 2 is to plot a histogram of the spike produced 
% in the 3 seconds following the first 10 flashes. 
% Then plot a histogram of the of the spikes produced by flahses between 
% seconds 240-360.  Repeat for flashes between 1600-1720 s.

%% STA Analysis

% In this section you will get introduced to analyzing STAs.  Let's begin
% by clearing the workspace

clear

% Let's now load an new data file with STAs.

datarun = load_data('/Volumes/Berlioz/Analysis/2014-09-23-0/data005/data005');
datarun = load_neurons(datarun);
%%datarun = load_data('/Volumes/lab/Experiments/Array/Analysis/2013-11-27-0/data001/data001');


% now for something new...
% load the STAs for all cells
datarun = load_sta(datarun, 'load_sta', 'all');

% you should also check out the help for load_sta to understand how it
% works.

help load_sta

% Try setting the parameter "verbose" to true. See what happens.  If you do
% not know how to changes these parameters, have a look at the inputParser
% helpf for matlab. 
% http://www.mathworks.com/help/matlab/ref/inputparser-class.html

% now load the params file.
datarun = load_params(datarun);

% again try setting the verbose optional input to true!  What does it tell
% you?

datarun = load_params(datarun, 'verbose', true);

% OK, now you have loaded the neurons file, the sta file and the params
% file into the datarun structure.  What fields does load_params add to the
% datarun?  What fields does load_sta add to the datarun?

% Check out the following
STA_size = size(datarun.stas.stas{1})

% The size is 60x60x1x32!  This is a big array!  That is why it took a lot
% of time to load the STAs for all the neurons.  The STA is 60 stixels
% wide, 60 high, 1 color, and 32 frames long.  The last frame is the frame
% that corresponded to the time of the spike. 

% Check out
RF_size = size(datarun.stas.rfs{1})

% Notice that this is empty and has a size of 0. 

datarun.stas

% notice that this structure contains many fields: stas, rfs, marks,
% rf_com (centers of mass), time_courses, etc.  Notice that all of these
% mostly empty.  That is becasue the receptive field (RF) and time courses
% (TCs)  rf_coms, and marks need to be calculated.  There are several ways
% to calculate these, but the simplist is to use the following function.

datarun = get_sta_summaries(datarun, 'all');

% Exercise: FIGURE OUT WHAT GET_STA_SUMMARIES IS DOING. Write a
% description.  Remember, don't edit this or any of these functions unless
% you make a copy and edit the copy only!  We all share these functions and
% changing them can impact or break someone else's analysis.

% you should run
help get_sta_summaries
% and read this carefully!

% Now look to see what is in datarun.stas.rfs, datarun.stas.marks,
% datarun.stas.time_courses, etc.

% Note that every cell has a cell ID assoicated with it.  You can see these
% in:
datarun.cell_ids

% the ID of the first cell is 17.  IN MOST FUNCTIONS, CELLS ARE REFERENCED
% BY THEIR ID #!  For example.

plot_rf(datarun,17)

%or

plot_rf(datarun, datarun.cell_ids(1))

% Notice this doesn't work!

plot_rf(datarun,1)

% to plot the receptive field (RF) of the first cell you need to refer to
% it by its cell_ID, not by its cell_index (ie 17, not 1).  This is
% important to remember and keep track of when you need to use an ID #
% versus an index.

% In Vision, I cell typed some of these RGCs.  You can see these cell types
% in

datarun.cell_types

% you can also run...

info(datarun)

% to see the different cell types, and their cell_ids.

% Let's plot the receptive field mosaics fpr ON type1 cells

% make vision fits the default fits to use
datarun = get_sta_fits_from_vision(datarun);

% plot the mosaic of recepive fields
plot_rf_summaries(datarun, 'ON type1', 'plot_fits', true)

% You should help and study

help plot_rf_summaries % This function can do a lot!

% Exercises: Add cel ID labels.  Make the labels red.  Plot the fits in blue
% Plot the mosaic of receptive fields for the OFF type1 cells.
% Calculate the mean RF diameter of the ON type1 cells and compare this to
% the OFF type1 and ON type3.

% Now let's look at the time courses.

help plot_time_courses

plot_time_courses(datarun, 'ON type1', 'bw', true)

% you can also do
plot_time_courses(datarun, 'ON type1', 'all', true, 'bw', true)
% to see them superimposed.

% notice that one of the time courses is pretty noisy and doesn't look like
% the other ones.  Let's see if we can fix this.

marks_params.thresh = 4;
datarun = get_sta_summaries(datarun, 'all', 'marks_params', marks_params);

plot_time_courses(datarun, 'ON type1', 'all', true, 'bw', true)

% plot noisy time course went away!  Exercise: What did I just do to clean
% up the time course?  Explain.























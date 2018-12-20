%% Create epochs, epochGroups, and populations.

% ------------------------------------------------------------------------------
% Import data into epochs.

% Create 50 empty epoch objects and store them in an array called 'epochs'.
epochs(50) = Epoch();

% Add response data to epoch. We'll use random data as a response for this 
% example. This is where you would import your data from file.
for i = 1:length(epochs)
    epochs(i).response = rand(100, 1);
end

% Now we have 50 epochs stored in an array called 'epochs'.

% ------------------------------------------------------------------------------
% Form epochGroups (usually cells in your case).

% Create 3 empty epoch group objects and store them in an array called
% 'epochGroups'.
epochGroups(3) = EpochGroup();

% Add appropriate epochs to cell 1.
epochGroups(1).epochs = epochs(1:10);

% Add appropriate epochs to cell 2.
epochGroups(2).epochs = epochs(11:20);

% Add appropriate epochs to cell 3.
epochGroups(3).epochs = epochs(21:50);

% Now we have 3 cells (epoch groups), containing our 50 epochs, stored in
% an array called 'epochGroups'.

% ------------------------------------------------------------------------------
% Form populations.

% Create 2 empty population objects and store them in an array called
% 'populations'.
populations = Population();

% Population 1
populations(1).epochGroups = epochGroups(1:2);

% Population 2
populations(2).epochGroups = epochGroups(1:3);

% Now we have 2 populations, containing different subsets of our 3 cells, stored
% in an array called 'populations'.


%% Perform analysis.

% Perform epoch level analysis.
for i = 1:length(epochs)   
    calcMean(epochs(i));
    %calcNumSpikes(epochs(i));
    %etc...
end

% Perform epoch group level analysis.
for i = 1:length(epochGroups)
    calcMean(epochGroups(i));
end

% Perform population level analysis.
for i = 1:length(populations)
    calcMean(populations(i));
end


%% Extract results.

% We'll just display them but it's just as easy to put them into a struct
% to transfer them to Igor.

for i = 1:length(epochGroups)
    disp(['Cell ' num2str(i) ' Mean: ' num2str(epochGroups(i).results.mean)]);
end

for i = 1:length(populations)
    disp(['Population ' num2str(i) ' Mean: ' num2str(populations(i).results.mean)]);
end


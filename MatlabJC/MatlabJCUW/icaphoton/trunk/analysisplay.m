% Adding Analysis to Playtime:
%
% Create an M-file that does the analysis desired.  This can be a separate m-file,
% or it can be a new case within Playtime.  The advantage of creating a separate m-file 
% is that it is a simple way to distribute new plotting commands to different machines,                      
% the disadvantage is that the easiest method of passing the variables back and forth is, 
% as of yet, still in question.  Regardless, the Playtime code needs to be changed to allow      
% access to the new command, so it may be easiest to save the code as an m-file for           
% distribution, and than copy the code into Playtime as a new case on various machines as needed.  
% The epoch data itself is stored within Playtime as a cell array called groupepochs.  
% The first row of cells is the data, and the first column of the second row is a list of 
% the epoch numbers.  Here is an example of how the epoch data is accessed to create a 
% histogram of the maximums
%
%		case 'maximum'			
%			want = size(groupepochs,2);
%			condense = cat(1,groupepochs{1,1:want});
%			condens = abs(condense)';
%			histoplot = max(condens);
%			playtime histogram
%	
% The cat command takes the cells in groupepochs and converts them to a single 
% matrix.  The abs command takes the absolute value.  If the completed analysis
% is to be in the form of a histogram, assign the analysis to the variable histoplot
% and than call playtime histogram.
%
% The next thing is to add a call to the menu for your new analysis.  At the 
% beginning of Playtime (somewhere around line 120) there are two blocks of code 
% that start:
%
%		h1 = uimenu('Parent',h0, ...
%
% The first of these is reserved for analysis that will be plotted on the same 
% plot as the data.  If there is interest in learning how to do this type of 
% analysis plotting, I will write a help file.
%
% The second chunk of code is to add a new graph to the desktop, specifically a 
% histogram. To add a histogram plot option to the histogram menu, add a line to 
% this block of code with the same syntax as the ones already there, ie.
%
%		uimenu(h1,'Label','Maximum','Callback','playtime maximum');
% 
% Where you will enter your own label according to how you want the analysis called,
% and the callback is the name of your m-file (one word) or the case added to 
% playtime (playtime newcase)
%  
% It is, of course possible to do other plots, for instance pie plots or 3-d plots
% 
% and it would be preferable to set up a new menu for these.  This is a very straight-forward
% process of copying the code of one of the other menus.
% 

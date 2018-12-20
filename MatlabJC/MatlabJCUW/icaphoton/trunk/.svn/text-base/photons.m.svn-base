% Icaphoton
% Version 2	OS X	15 April 2002
% 
% Icaphoton (I-see-a-photon) is a collection of MATLAB functions (m-files) used 
% to retrieve data from the original data files created in Igor, put it into a 
% matlab structure (CellINfo), and analyze it.  Uses 'ITC' mex files from Aquirino 
% Extensions to retrieve data.  If you have problems with crashing,
% double-check that you have the latest mex files (found here on garfield:
% Users/lab/Public/Local/Programming/Acquirino\ version\ 6/OSX-mex-files-version7)
% 
% other help files:
% help analysisplay	   Explains how to add additional analysis tools to Playtime
% help StructCellInfo    Explains the matlab data structure 
% 
% Each of the following files also contains help information, type help filename
% 
% Main GUIs (Graphical User Interface):
% icaphoton         GUI. First interface.  Used to select data file.
% epochs            GUI. Second interface.  Reached by going through icaphoton.
%                       Used to select data set by choosing all data that
%                       is exactly the same by some user-defined criteria.
%                       Can also go to familysearch from here using button
%                       'Find Families'. See help file
% playtime          GUI. Third interface.  Can fine-tune data selection and
%                       do some preliminary analysis.
% familysearch      GUI. Select data in families sets.  Families are groups
%                       of data that increase in time in some feature, most
%                       often by amp, but can be any feature.  Must have
%                       more than two epochs per family for the family to
%                       be selected. 
% 
% Non-GUI functions related to CellInfo structure:	
% LoadCIData        function to load epoch data into the structure CellInfo 
% LoadSCIData       version of LoadCIData for data files with more than one
%                       segment
% LoadCIStim        function to load the stimulus data into the structure
%                       CellInfo
% RmCIData          function to remove epoch data (either type) from CellInfo
% FindFamPara       function to find the value of a parameter in a family
%                       condition structure in CellInfo. Looks in fields
%                       FamilyCues and FamilyStep.
% FindSearchPara    function to find the value of a parameter in a
%                       condition (can be either family or epoch) structure
%                       in CellInfo.  Looks in field SearchCrit
% GetCellData       takes an array of cell numbers (epochnumber + 1, since
%                       there is no zero cell in matlab), and creates a
%                       matrix of the data. Since it will be in matrix
%                       form, the epoch data for these epochs must all be
%                       the same length.
% makenewfactor     script created to convert an old-style CellInfo in
%                       which the fields OutputScaleFactor and
%                       NDFConfiguration contain strings to ones that have
%                       numbers.  Will go through every file in the CURRENT 
%                       directory and check to see if the mat files are
%                       ones that contain the old-style fields.  Should not
%                       do anything to either newer CellInfo files or files
%                       that are not mat or do not contain CellInfo
% 
% changeFilePath        script to run in a directory that has mat files with
%                       the path names in os 9 format to change the paths
%                       to os X format.  Will not bother other files.
% 
% Calibration functions: ('help calibration' for more information)
% 
% saveAllSpectra    Prior to calibrating, load the appropriate spectra
%                       (both animal and cell type) by using this script.
%                       A file called allSpectra will be created.  This
%                       file contains the spectra for both the LEDs and the
%                       cells. Make sure the allSpectra.mat file it creates
%                       is on your path, becasue the other calibration
%                       routines will need it. 
% LEDCalibrate      Calibrates a file.  Interactive, shows plots of spectra.  Puts the
%                       calibrated stimulus amp into the sub-field
%                       StimulusAmp in the UserInfo field of the
%                       appropriate Condition.  For FamilyConditions if the
%                       FamilyStep field is StmAmp, it also makes a matrix
%                       the size of all epochs in the condition with the
%                       calibrated StmAmp for each epoch. This is also in
%                       the field UserInfo, but in the subfield FamilyStep.  
% AutoLEDCalibrate	Does the same as LEDCalibrate, but not interactive, allows you 
%                       to enter LED Setting and Spot size as inputs
%                       instead (assumes same for whole file), and
%                       automatically goes through all family and epoch
%                       conditions.  Does not plot.  
% batchcalibrate	Allows user to calibrate an entire CURRENT directory automatically.  
%                       Assumes all files had same settings for LED and spot size.  Must
%                       edit the script itself to set these.  Should do nothing to a file
%                       that does not have a CellInfo or is not a mat file.
% 
% Family functions:
% 
% GUI interface:
% familyplay		Plots families for analysis.  Returns the x and y coordinates of 
%                       each plot to the command window, where it can be saved to a file, 
%                       either for further analysis or to export to igor or some other 
%                       graphing program.  Various choices on types of plots
% 				
% Non-GUI functions:
% FamilyMax       	Returns a max for all the maxes of epochs in this 
%                       familycondition. Useful for normalizing.
% GetFamilies		Gives back a list of epochnumbers and the families to which 	
%                       they belong, that correspond to the steps and cue given.  
% ParseFamilies   	Gives back data and epochnumbers parsed into families (separate 
%                       cells for each family, matrix in each cell where each row is a 
%                       different epoch), given a list of epcohnumbers and the families
%                       to which they belong (ie. output of GetFamilies).
% CueMax			Input CellInfo, number of the family you are interested in and the
%                       cue (row) number. Cues are the parameters that can change between 
%                       families - ie. maybe you have some families with a red led and 
%                       some with a green led and if you want the maxes for all the 
%                       families based on the leds you must enter the corresponding 
%                       number for the cues, this number is the row number in the matrix 
%                       FamilyCues. The function will then return the max for each of 
%                       these cues. Userful for normalizing.
% datamaxorarea		Given a data set from familyplay (xdata and ydata from the 
%                       command window, data sits in cells within these variables) finds 
%                       the max and area for all epochs in each cell (one max and area per 
%                       cell).
% familymodel     	Used with nlinfit to model the data plotted R/Rmax vs. Flash 
%                       Strength.  See help nlinfit.
% famindex        	Given epochnumbers finds the index numbers in the FamilyCondition 
%                       user inputs that corresponds to these numbers.
% familymodel 		function that is used to fit the Sensitivity vs. log of flash
%                       intensity.
% fitroutine		takes xvalue and yvalue loaded into workspace from familyplay
%                       Should be a whole group of families, not averaged.  Will run 
%                       nlinfit with fitmodel for each family, and plot the average, with
%                       a fit.  Eventually will plot error bars using standard error, and 
%                       run a loop to get a better fit.  Work in progress.
% 

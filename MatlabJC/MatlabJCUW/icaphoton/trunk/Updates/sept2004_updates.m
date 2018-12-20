% On some machines, the ITCFinishAnalysis in playtime before the user is sent
% back to the epochs interface causes matlab to crash.  In the latest version of 
% playtime, I have commented out this function call until the cause of this can
% be determined.

% It isn't terribly clean to not use the finish analysis, but it causes matlab crashes on some 
% machines.  The analysis is finished in every other interface, including this one when we 
% leave for good, so this shouldn't be a big deal.  Why it only does it here is truly a mystery.

% Another change is returning the functionality to loading CellInfo to icaphoton.
% Ie., now you can go back to an old CellInfo and use ICAphoton to add more conditions.

% RmCIData had an error in it, which was corrected.

[s,r] = system('df | grep /Volumes/Groups | grep ^afp | awk ''{print $6}''') ; % this line finds the most recent "/Volumes/Group"
RootDir = [r(1:end-1),'/Lab Staff/Calibration/'] ;
cd(RootDir)
calibrationGUI

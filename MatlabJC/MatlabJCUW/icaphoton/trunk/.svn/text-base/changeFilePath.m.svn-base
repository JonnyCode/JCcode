% Run script in the directory that contains the mat (CellInfo) files that are still in
% os 9 format.  Assumes all mat files correspond to data files that have
% been converted to OS X, and are all in the same directory (which is not
% necessarily the current directory).  Assumes name of converted file is
% the same as the name of the original data file.  Does not bother files
% that are not mat files or that have paths that are already os x format.
% Created MKMK 12/04

newpath = input('What is the path to the CONVERTED (for OS X) data files? ','s');
if newpath(end) ~='/'
    newpath = [newpath '/'];
end
dircont = dir('*.mat');
for i = 1:length(dircont)
    filename = dircont(i).name
    load(filename)
    s = who;
    if any(strcmp(s,'CellInfo'))
        fp = CellInfo.CellFile;
        if ~isempty(findstr(':',fp))
            if ~strcmp(version,'5.2.1.1421')
                colplace = findstr(':',fp);
                filename = fp(colplace(end)+1:end);
                CellInfo.CellFile = [newpath filename]
                save(filename,'CellInfo')
            else
                disp('This mfile is not designed for OS 9')
            end
        else
            disp('This matfile does not contain an os 9 path')
        end
        clear CellInfo
    else
        disp('Did not change file')
        clear s
    end
end
clear



function waitForFileExist(dirname,fname)
%dirname
%fname
curDir = pwd;
cd(dirname)
while ~exist(fname,'file')
%    1
    %do nothing
    WaitSecs(0.05); %50 ms
end
%delete(fname);

cd(curDir);





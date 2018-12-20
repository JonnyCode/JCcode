%function out=readIgor(filename,type)

%filename - a string

%type - a string that designates the type of number stored in the wave.

%       the only ones supported are 'float32' and 'int16' but adding others

%       is easy. just make sure you take off the correct number of elements

%       at the end and use fread in the correct format

%

%you have to change the fread statement if you are reading files that were

%created on a pc

%

%this may not be fancy efficient c-code but its pretty fast and convenient!

%-Andrew Gartland

function out=readIgor(fn,type),

fid=fopen(fn,'rb');



%skips the 126 bytes of header

fseek(fid,126,'bof');



%reads the remainder of the file

%'ieee-be' denotes the byte order as big-endian (unix and mac)

%change to 'ieee-le' for files that were created on a pc with little-endian

out=fread(fid,inf,type,'ieee-be');

fclose(fid);



%take off last 16 bytes!

if(strcmp(type,'int16')),

    out=out(1:length(out)-8);

elseif(strcmp(type,'float32')),

    out=out(1:length(out)-4);

end

end


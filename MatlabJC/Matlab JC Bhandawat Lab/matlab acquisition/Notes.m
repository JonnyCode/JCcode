function FullNote = Notes(filePath)

% for writing notes during an experiment
% JC 4/16/12

Note = input('Notes:','s') ;

NoteTime = datestr(now) ;

FullNote = [NoteTime,' Note: ',Note] ;

fid = fopen(filePath,'a');
fprintf(fid,'%s \r\n',FullNote);
fclose(fid);






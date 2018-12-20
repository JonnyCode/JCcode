function s = cellArraysToMatricesInStruct(s)
fnames = fieldnames(s);

for i=1:length(fnames)
   if iscell(s.(fnames{i}))
       s.(fnames{i}) = cell2mat(s.(fnames{i}));
   end    
end
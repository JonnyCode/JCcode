function DataStruct = hdf52Struct(FilePath) 

%This function takes an hdf5 file and makes a structure.  Made this
%function with hdf5 files created from igor experiments in mind.

%JC 1/2/12

%FilePath = string with file path/filename
%FilePath = C:\Cafaro Data\111230Data\111230#1.h5 ; example file path

DataStruct = struct('FilePath',FilePath) ; % set structure and file path

Info = hdf5info(FilePath) ; % get file info

for a=1:length(Info.GroupHierarchy.Datasets) ; % for every dataset in the hdf5
    TempData = hdf5read(Info.GroupHierarchy.Datasets(a)) ; % grab the data
    TempFieldName = Info.GroupHierarchy.Datasets(a).Name ; % filed name
    FieldName = TempFieldName(2:end) ; % get rid of first character since it seems to always be a '/'
    DataStruct = setfield(DataStruct,FieldName,TempData) ; % set field in structure
end



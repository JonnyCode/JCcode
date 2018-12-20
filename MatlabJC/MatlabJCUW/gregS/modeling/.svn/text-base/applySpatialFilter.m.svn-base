function M = applySpatialFilter(M)
%M is th FFModel
%applies spatial filter specified in layer 1 
%to stimulus then sums across space and fills in iputs to
%layer 1 subunits

T = TwoD2OneDTransform;

%to try for speedup: sparse matrices
nUnits = length(M.layers(1).subunits);

S = M.stimulus; %2D in space, 1D in time, time is z
[Lx, Ly, Ltime] = size(S);
disp('Applying spatial pre-filters...')
for j=1:nUnits %for each subunit
    disp(['Unit ' num2str(j) ' of ' num2str(nUnits)]);
    V = zeros(1,Ltime);
    for i=1:Ltime %for each frame
       V(i) = T.invert(squeeze(S(:,:,i))) * T.invert(M.layers(1).subunits(j).spatialPreFilter)';       
    end
    M.layers(1).subunits(j).input = M.layers(1).transform.invert(V);    
end

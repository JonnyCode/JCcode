function M = generateTexture2(l,w,sigma,C)
% modified previous function to make rectangle instead of square
% JC 6/25/10
%C is contrast

M = rand(w,l);
if sigma>0
    win = fspecial('gaussian',w,sigma);
    M = imfilter(M,win,'replicate');
    M = M./max(M(:));    
else
    %do nothing
end

M = makeUniformDist(M,C);

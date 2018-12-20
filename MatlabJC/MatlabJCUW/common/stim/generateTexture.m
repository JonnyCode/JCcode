function M = generateTexture(S,sigma,C)
%C is contrast
winL = 200;
M = rand(S,S);
if sigma>0
    win = fspecial('gaussian',winL,sigma);
    win = win ./ sum(win(:));
    M = imfilter(M,win,'replicate');
    M = M./max(M(:));    
else
    %do nothing
end

M = makeUniformDist(M,C);

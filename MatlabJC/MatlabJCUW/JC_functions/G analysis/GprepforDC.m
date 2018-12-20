% this function will take recorded conductances and prepare them for
% dynamic clamp

for a=1:size(G,1) ;
    GmeanSub(a,:) = G(a,:)-mean(G(a,:));
end

minG = min(GmeanSub(1:end)) ;
Goffset = GmeanSub-minG ;

maxG = max(Goffset(1:end)) ;
GforDC = Goffset/maxG ;


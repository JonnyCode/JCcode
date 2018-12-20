% for response to reviewers question about spike detection in PNs

% 11/7/2015 JC

figure 
n=0 ;

for a = 1:NumBackgrounds ; % for each background
    for b = 1:NumConcentrations ; % for concentration
        for c = 1:NumTrials(a,b) ;
            n=0 ;
            pnt = spikePnt{a}{b}{c} ;
            for d = 1:length(pnt) ;
               
                if pnt(d)>500 && pnt(d)<length(vData{a}{b}(c,:))-500 ;
                    n=n+1 ;
                    spikeData{a}{b}{c}(n,:) = vData{a}{b}(c,pnt(d)-500:pnt(d)+500)-vData{a}{b}(c,pnt(d)-40) ;
                end
            end
        end
    end
end
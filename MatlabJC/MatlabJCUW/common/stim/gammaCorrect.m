function monitorLevel = gammaCorrect(I, gammaTable)
%I is from 0 to 1
gammaTable = gammaTable - min(gammaTable); %zero offset
g_norm = gammaTable./max(gammaTable);
monitorLevel = zeros(size(I,1),size(I,2));

if max(I(1:end))>1 ;
    I(I>1)=1 ;
    disp('I > 1')
    
elseif min(I(1:end))<0 ;
    I(I<0)=0 ;
    disp('I < 0')
end

for i=1:numel(I)
    monitorLevel(i) = find(g_norm>=I(i),1)-1;%0 to 255
end

 

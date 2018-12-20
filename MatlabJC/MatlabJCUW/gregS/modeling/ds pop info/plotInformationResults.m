function [] = plotInformationResults(T,pairs,testAngles,I2,I2_corr)
Npairs = size(pairs,1);
phaseOffset = 360./size(T,1);

%get all 'unique shift pairs'

D = rem(phaseOffset.*abs(pairs(:,1) - pairs(:,2)),180);
D(D==0) = 180;
[D_unique, Ind] = unique(D);

pairs = pairs(Ind,:);
Npairs = size(pairs,1);
I2 = I2(Ind,:);
I2_corr = I2_corr(Ind,:);

X = 0:359;
for i=1:Npairs
    subplot(Npairs,2,(i-1)*2+1);
    plot(X,T(pairs(i,1),:),'b')
    hold on;
    plot(X,T(pairs(i,2),:),'g')
    %plot(X,T(pairs(i,1),:)+T(pairs(i,2),:),'k')
    hold off;
    subplot(Npairs,2,i*2);    
    plot(testAngles,I2(i,:),'kx-'); hold on; plot(testAngles,I2_corr(i,:),'rx-');
    axis([0 360 0 1]);
    title([num2str(D_unique(i)) ' degrees apart']);
    hold off;
end
hold off;
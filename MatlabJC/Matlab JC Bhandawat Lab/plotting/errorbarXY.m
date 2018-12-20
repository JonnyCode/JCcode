function  errorbarXY(Xmean,Ymean,Xerror,Yerror,pointString,colorString,figNum) 

% vector with points and errors
% JC 7/8/13

figure(figNum)
for a = 1:length(Xmean) ; % for each point on the graph
    plot([Xmean(a),Xmean(a)],[Ymean(a)-Yerror(a),Ymean(a)+Yerror(a)],'-','color',colorString,'LineWidth',1)
    hold on
    plot([Xmean(a)-Xerror(a),Xmean(a)+Xerror(a)],[Ymean(a),Ymean(a)],'-','color',colorString,'LineWidth',1)
    plot(Xmean(a),Ymean(a),pointString,'color',colorString)
end

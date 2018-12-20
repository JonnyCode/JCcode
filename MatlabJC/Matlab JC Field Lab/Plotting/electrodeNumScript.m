% Xiaoyang's code to plot electrode numbers

rotation = -90;
position = datarun.ei.position;

% rotates
for i = 1:length(position)
    [t1 t2]=cart2pol(position(i,1),position(i,2));
    t1=t1+(rotation/(180/pi));
    [position(i,1),position(i,2)]=pol2cart(t1, t2);
end

% plots
for ee=1:size(position,1)
    labeltext = num2str(ee);
    text(position(ee,1),position(ee,2),labeltext,'Color', 'k','FontSize',10,...
        'HorizontalAlignment','Center','VerticalAlignment','Bottom');
    hold on
end

xlim([-400 400])
ylim([-400 400])
axis equal

% some other code that does not require an ei

ai = load_array_info(1530) ; % for 519 30 um spacing array

e=[1,23,100] ; % electrodes
plot_electrodes(ai.positions,{e},'alpha', false) ; % plot array with electrodes e colored in











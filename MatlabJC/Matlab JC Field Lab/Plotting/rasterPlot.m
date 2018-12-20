%spt is a cell array where each element is a vector of spike times
%dur is a two element vector which allows you to specify a group of spike that should be plotted as a different color. right now its turned off. AJG, 10/27/08
function rasterPlot(spt,dur)
for i=1:length(spt)
    for j=1:length(spt{i})
        if nargin==2 && spt{i}(j)>dur(1) && spt{i}(j)<dur(2)
            col='k';
        else
            col='r';
        end
        plot([spt{i}(j) spt{i}(j)],[i i+1],col);
        hold on;
    end
end
hold off;
a=axis;
axis([a(1:2) 1 length(spt)+1]);
end
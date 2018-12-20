function [maxes,areas] = datamaxorarea(ydata,xtime)
% Given a data set from familyplay (xdata and ydata, sits in cells), finds the max
% and area for all epochs in each cell (one max and area per cell).
maxflag = 0;
ind = 0;
allmaxes = zeros(length(ydata),1);
allmaxarea = zeros(length(ydata),1);
for i=1:length(ydata)
  testmax = zeros(length(ydata{i}),1);
  testarea = zeros(length(ydata{i}),1);
  for j = 1:length(ydata{i})
	 test = ydata{i}{j};
	 testmax(j) = max(max(test));
	 for k = 1:size(test,1)
		newtempdata = test(k,:);
		newtempdata(1) = newtempdata(1)/2;
		newtempdata(end) = newtempdata(end)/2;
		I = xtime*sum(newtempdata);
		I = abs(I);
	 end
	 testarea(j) = I;
  end
  allmaxes(i) = max(testmax);
  allmaxarea(i) = max(testarea);
end
maxes = allmaxes;
areas = allmaxarea;

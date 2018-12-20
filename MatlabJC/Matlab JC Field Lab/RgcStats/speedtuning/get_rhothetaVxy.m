function[rho, theta, num_t, num_c] = get_rhotheta(NumSpikesCell2, StimComb2, datarun)

%Function extracts spike numbers and angles of each direction for each temporal period

%Inputs: NumSpikesCell2: Total avergae spike numbers for each trial 

        %StimComb2: List of all combinations of S & T periods & directions
        
%Outputs: rho: spike numbers of each direction for each temporal period

          %theta: angles of each direction for each temporal period
          
          %num: List of all temporal periods
          
          %Sneha Ravi 
  %Last revision: 12-18-2012

num_t = unique(cell2mat(StimComb2(:,2))); %Number of temporal periods
% num_c = unique(cell2mat(StimComb2(:,4)), 'rows'); %Number of contrast/bar_color
% num_t = datarun.stimulus.params.DELTA;
num_c = datarun.stimulus.params.RGB;
rho = cell(length(num_t), length(num_c));
theta= cell(length(num_t), length(num_c));
for i = 1:length(num_t)
    for j = 1:length(num_c)
        [r,t] = deal([]);
        int = ismember(cell2mat(StimComb2(:,2)), num_t(i))';
        inc = ismember(cell2mat(StimComb2(:,4)), num_c{j}, 'rows')';
        r = NumSpikesCell2(:,int&inc);
        t = cell2mat((StimComb2(int&inc,3)'))*pi/180;
        t = repmat(t,size(r, 1), 1);
        rho{i,j} = r;
        theta{i,j} = t;
    end
end
end


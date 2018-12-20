function[C] = eiplots(datarun, cell_ids, varargin)
%Function plots electrophysiological images in different ways with different parameters
%Input: Datarun structure with eis
%       Cell IDs of target cells

% 1) 'plotei': Plots EI with different parameters ('cutoff', 'alpha', 'scale', 'color_palette')
% 2) 'exclude_axon': Plots EI without the axon, with a 'threshold' value
% 3) 'plot_sourceneighb': Plots source elecrode and neighboring electrodes
% 4) 'fit_Gaussian': Fits EI with 2-D Gaussian
%                     Gaussian can be fit in 3 ways ('fit_type'):
%                     'gmdist': uses matlab function gmdistribution to construct Gaussian from given mean and covariance
%                     'fitgm': uses matlab function gmdistribution.fit to fit Gaussian - no initial conditions
%                     'fitgmic' uses matlab function gmdistribution.fit to fit Gaussian - with initial conditions
%                     Gaussian can be visualized in 3 ways ('contour_type'): 
%                    'simplecontour': uses matlab function contour, can adjust number of contour lines
%                    'ezcontour': uses matlab function ezcontour
%                    'sdcontour': plots gaussian fit n number of std deviations from mean

p = inputParser;
p.addParamValue('color_palette', []); % RGB color values for each cell's EI
p.addParamValue('cutoff', 0.03); % threshold for excluding electrodes with Voltage below cutoff
p.addParamValue('alpha', 0.5); % Intensity of color for each electrode
p.addParamValue('scale', 1); % Scales electrodes that are not the source electrode 
p.addParamValue('pause', false); % Allows plot to pause for each cell
p.addParamValue('save', false); % Saves image as a pdf to directory pathname as filename
p.addParamValue('pathname', '');
p.addParamValue('filename', '');
p.addParamValue('exclude_axon', false); % If you want to exclude axon from the EI
p.addParamValue('threshold', 3); % Cutoff threshold for axons
p.addParamValue('plotei', false); % Plots EI
p.addParamValue('plot_sourceneighb', false); % Plots source electrode and neighbors
p.addParamValue('neighbor_edgecolor', 'b');
p.addParamValue('neighbor_facecolor', 'c');
p.addParamValue('neighbor_size', 6);
p.addParamValue('fit_Gaussian', false); % Plots Gaussian fit to EI
p.addParamValue('fit_type', ''); % Whether to fit with initial conditions
p.addParamValue('contour_type', ''); %Way to plot the contour
p.addParamValue('contourgrid', 60); %Number of grids to use for ezcontour 
p.addParamValue('numcontours', 1); % Number of contours to display for simplecontour
p.addParamValue('stddev', 1); %Number of standard deviations for Gaussian fit

p.parse(varargin{:});
params = p.Results;

temp_indices = get_cell_indices(datarun, cell_ids);

%% Axon Included, EI with many parameters

for cc = 1:length(temp_indices) 
    temp_ei = datarun.ei.eis{temp_indices(cc)};
    if params.exclude_axon
        res =  ei_get_axon(datarun, cell_ids(cc), 'threshold', params.threshold); %Get electrode numbers for axon
        temp_ei(res,:) = 0;
    end
    if params.plotei
        plot_ei_(temp_ei, datarun.ei.position, 0, 'cutoff', params.cutoff, 'elec_colors', repmat(params.color_palette(cc,:), 512, 1), 'alpha', params.alpha, 'scale', params.scale);
        hold on;      
    end
    if params.plot_sourceneighb
        ei = temp_ei/max(abs(temp_ei(:)));
        %ei(abs(ei)<1) = 0; %Set all electrodes except source electrode with maximum voltage to zero
        [r c] = find(abs(ei)==1); %Find source electrode
        %pos  = datarun.ei.position(get_ei_neighbors(r,512),:); %Get neighboring electrode positions
        pos = get_ei_neighbors(r,512); %Get neighboring electrode positions
        ei(setdiff( 1:1:512, pos),:) = 0; %Only have voltage values of electrode and neighbors
        %plot(pos(2:end,1), pos(2:end,2), 'o','MarkerEdgeColor',params.neighbor_edgecolor,'MarkerFaceColor',params.neighbor_facecolor, 'MarkerSize',params.neighbor_size); %neighbors are equal size
        hold on;
        plot_ei_(ei, datarun.ei.position, 0, 'cutoff', params.cutoff, 'elec_colors', repmat(params.color_palette(cc,:), 512, 1), 'alpha', params.alpha, 'scale', params.scale);
    end
    if params.fit_Gaussian
        [x y] = meshgrid(datarun.ei.array_bounds_x(1,1):1:datarun.ei.array_bounds_x(1,2), datarun.ei.array_bounds_y(1,1):1:datarun.ei.array_bounds_y(1,2));
        ei = temp_ei/max(abs(temp_ei(:)));
        ei(abs(ei)<params.cutoff) = 0; %Set all electrodes below cutoff to zero
        ei_frame = get_ei_max_frame(ei, 1);
        elecofei = find(ei_frame);
        ei(abs(ei)<1) = 0;
        [r2 c2] = find(ei);
        elec_cent = r2;
        allx = datarun.ei.position(elecofei,1); %get all x and y positions of electrodes above cutoff
        ally = datarun.ei.position(elecofei,2);
        X = [];
        X(:,1) = allx;
        X(:,2) = ally;
        mu = [datarun.ei.position(elec_cent,1), datarun.ei.position(elec_cent,2)]; %Source electrode is the mean
        if(strcmp(params.fit_type, 'gmdist'))
            obj = gmdistribution(mu,cov(X));
        elseif(strcmp(params.fit_type, 'fitgm'))
             obj = gmdistribution.fit(X,1, 'Replicates',1,'CovType', 'full'); % If you want to fit the EI - no initial conditions
        elseif(strcmp(params.fit_type, 'fitgmic')) %Fit EI with initial conditions
            %S.mu = mu;
            %S.mu = mean(X);
            mass = abs(ei_frame(ei_frame~= 0));
            pts = centroid(X, mass);
            S.mu = pts;
            S.Sigma = cov(X);
            obj = gmdistribution.fit(X,1, 'Replicates', 1, 'CovType', 'full', 'Start', S);
        end  
        plot(datarun.ei.position(elec_cent,1), datarun.ei.position(elec_cent,2), '+r');
        hold on;
        %plot(pts(1), pts(2), '+b');
        if (strcmp(params.contour_type, 'ezcontour'))
            ezcontour(@(x,y)pdf(obj,[x y]),[datarun.ei.array_bounds_x(1,1) datarun.ei.array_bounds_x(1,2)],[datarun.ei.array_bounds_y(1,1) datarun.ei.array_bounds_y(1,2)], params.contourgrid);
            hold on;
        elseif (strcmp(params.contour_type, 'simplecontour'))
            for a = 1:size(y,2)
            YY(:,a) = mvnpdf([x(:,a) y(:,a)], obj.mu, obj.Sigma); %List of probabilities at each position in the 2-d gaussian
            end
            contour(x,y,YY, params.numcontours);
            hold on;
        elseif (strcmp(params.contour_type, 'sdcontour'))
            for a = 1:size(y,2)
                YY(:,a) = mvnpdf([x(:,a) y(:,a)], obj.mu, obj.Sigma); %List of probabilities at each position in the 2-d gaussian
            end
            if (round(obj.mu(1,1)+ params.stddev*std(allx)) < datarun.ei.array_bounds_x(1,2)) %Check if standard deviation limit is within bounds of the array for plotting
                [rx cx] = find(x==round(obj.mu(1,1)+ params.stddev*std(allx)));
            else
                [rx cx] = find(x==round(obj.mu(1,1)- params.stddev*std(allx)));
            end
            if (round(obj.mu(1,2)+ params.stddev*std(ally)) < datarun.ei.array_bounds_y(1,2))
                [ry cy] = find(y==round(obj.mu(1,2)+ params.stddev*std(ally)));
            else
                [ry cy] = find(y==round(obj.mu(1,2)- params.stddev*std(ally)));
            end
            v = [YY(ry(1,1),cx(1,1)), YY(ry(1,1),cx(1,1))];
            [C h] = contour(x,y,YY,  v, 'b', 'LineWidth', 2);
        end
    end
    if params.pause
        pause; %Pause plot for each cell
    end
end

axis off

if params.save
    save_figure_pdf(params.pathname, params.filename, gcf) %Save the figure
end

end
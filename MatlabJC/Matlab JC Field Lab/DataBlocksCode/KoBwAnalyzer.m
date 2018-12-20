function ForIgor = KoBwAnalyzer(DataBlock, DB, Params)

% this function will analyze binary white stimuli response from array data
% before and after cells are knocked out (KO)

% JC 11/11/14
DB=2; % TEMP

% parameters
mapEiFlag = DataBlock(DB).mapEiFlag ; % if datablocks should be mapped using map_ei

% load data
PreKo = load_data(DataBlock(DB).PreKo.DataPath) ;
PreKo = load_neurons(PreKo) ;

PostKo = load_data(DataBlock(DB).PostKo.DataPath) ;
PostKo = load_neurons(PostKo) ;

% how different are the preko and postko response fits
% spatial rf?
marks_params.thresh = 3 ;
PreKo = load_sta(PreKo) ;
PreKo = load_params(PreKo) ;
PreKo = get_sta_summaries(PreKo, 'all','marks_params', marks_params) ;
PreKo = get_sta_fits_from_vision(PreKo) ;

PostKo = load_sta(PostKo) ;
PostKo = load_params(PostKo) ;
PostKo = get_sta_summaries(PostKo, 'all','marks_params', marks_params) ;
PostKo = get_sta_fits_from_vision(PostKo) ;

% get list of cell types in PreKo
for a = 1:length(PreKo.cell_types) ;
    PreKo_celltypes{a} = PreKo.cell_types{a}.name ;
end

UniqueCellTypes = unique(PreKo_celltypes) ;
if isempty(UniqueCellTypes{1}) ;
    UniqueCellTypes = UniqueCellTypes(2:end) ;
end

% if using map ei cells
if mapEiFlag ;
    % load electrical images
    PreKo = load_ei(PreKo, 'all') ;
    PostKo = load_ei(PostKo, 'all') ;

    % map using electrical images
    [cell_list_map, failed_cells] = map_ei(PreKo, PostKo) ;

    % cells ids in PostKo for each UniqueCellType set in PreKo
    for a = 1:length(UniqueCellTypes) ;
        Tempi = get_cell_indices(PreKo, UniqueCellTypes{a}) ;
        post_ids{a} = cell2mat(cell_list_map(Tempi)) ;
    end
else
    for a = 1:length(UniqueCellTypes) ;
        post_ids{a} = intersect(PostKo.cell_ids, get_cell_ids(PreKo,UniqueCellTypes{a})) ;
    end
end

% plot pre and postKo rfs

for a = 1:length(UniqueCellTypes) ;
    figure(a) 
    subplot(1,2,1)
    plot_rf_summaries(PreKo, UniqueCellTypes{a}, 'plot_fits', true, 'foa',a,'fit_color','k') ; % plot mosiac
    
    hold on
    if ~isempty(post_ids{a}) ;
        plot_rf_summaries(PostKo, post_ids{a}, 'plot_fits', true, 'foa',a,'fit_color', 'r') ; % plot mosiac
    end
    title(UniqueCellTypes{a})
end

% temporal rf?
for a = 1:length(UniqueCellTypes) ;

    Tempi = get_cell_indices(PreKo, UniqueCellTypes{a}) ;

    figure(a) 
    subplot(1,2,2)
    for b=1:length(Tempi) ;
        plot(flipud(PreKo.stas.time_courses{Tempi(b)}),'k')
        hold on
    end

    
    if ~isempty(post_ids{a}) ;
        Tempi = get_cell_indices(PostKo, post_ids{a}) ;

        figure(a) 
        subplot(1,2,2)
        for b=1:length(Tempi) ;
            plot(flipud(PostKo.stas.time_courses{Tempi(b)}),'r')
            hold on
        end
    end
    
    title(UniqueCellTypes{a})
end


% how stable are the response fits for a single conditions?
% spatial rf?

% temporal rf?


% how well are the models at predicting the response?




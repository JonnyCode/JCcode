function raster_mb_all = combine_repeats(raster_mb)

raster_mb_all = cell(length(raster_mb), 1);
for cc = 1:length(raster_mb)
    if ~isempty(raster_mb{cc})
        for i = 1:size(raster_mb{cc}, 1)
            for j = 1:size(raster_mb{cc}, 2)
                for k = 1:size(raster_mb{cc}, 3)
                    for m = 1:size(raster_mb{cc}, 4)
                        raster_mb_all{cc}{i,j,k,m} = cell2mat(squeeze(raster_mb{cc}(i,j,k,m,:)));
                    end
                end
            end
        end
    end
end

end
function plot_shade(solution,modeColors,min_v,max_v)

    uniqueModes = unique(solution.mode);
    
    for i = 1:length(uniqueModes)
        mode = uniqueModes(i);
        color = modeColors(mode{:});
        idx_range = 1:(length(solution.Time));
    
        idx = idx_range(strcmp(solution.mode,mode));
        if idx(end) == idx_range(end)
            idx = idx(1:end-1);
        end
        
        for l = 1:length(idx)
            j = idx(l);
            fill([solution.Time(j)  solution.Time(j+1)  solution.Time(j+1) solution.Time(j)], ...
             [min_v min_v max_v max_v], ...
             color, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
        end
    end



end
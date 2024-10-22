function plot_shade_name(solution,min_v,max_v)
    pos_island = mean(solution.Time(strcmp(solution.mode,'Island')));
    text(pos_island, 1.0025*max_v, 'Island', 'Color', 'k', 'FontSize', 15,'HorizontalAlignment', 'center');
    

    ar_time_con = solution.Time(strcmp(solution.mode,'Grid connected'));
    [val, idx] = max(diff(ar_time_con));

    pos_connect = mean(ar_time_con(1:idx));
    text(pos_connect, 1.0025*max_v, 'Connected', 'Color', 'k', 'FontSize', 15,'HorizontalAlignment', 'center');
    
    pos_connect = mean(ar_time_con(idx+1:end));
    text(pos_connect, 1.0025*max_v, 'Connected', 'Color', 'k', 'FontSize', 15,'HorizontalAlignment', 'center');
   

    pos_resync = mean(solution.Time(strcmp(solution.mode,'Prepare for connection')));
    text(pos_resync, 1.0025*max_v, 'Resync', 'Color', 'k', 'FontSize', 15,'HorizontalAlignment', 'center');
    
    pos_prepare = mean(solution.Time(strcmp(solution.mode,'Prepare Island')));
    text(pos_prepare, 1.0025*max_v, 'Prepare', 'Color', 'k', 'FontSize', 15,'HorizontalAlignment', 'center');
end
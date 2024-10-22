
function f = make_plot_time2(solution,solution_connected,solution_island,Grid_para,constraints,range,modeColors)

set(groot,'defaultFigureVisible','on')

set(0, 'DefaultTextInterpreter', 'Latex')
set(0, 'DefaultLegendInterpreter', 'Latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'Latex')

colorblindSet = {[215,25,28]/255;[253,174,97]/255;[171,217,233]*.8/255;[44,123,182]*.8/255};


f = figure('Renderer', 'painters', 'Position', [10 10 850 400]);
    hold on 
        data_t_con= solution_connected.duration;
        data_t_isl= solution_island.duration;
        
        min_v = 0;
        max_v = 2;
    
        plot_shade(solution,modeColors,min_v,max_v)
        plot_shade_name(solution,min_v,0.95*max_v)
    
        p1 = plot(solution_connected.Time,data_t_con,'LineWidth',2,'Color',[0.2 0.2 0.2]);
        p2 = plot(solution_island.Time,data_t_isl,'LineWidth',2,'Color',colorblindSet{2});
         legend([p1 p2],'time connected','time island','Location','southwest','FontSize',15);
        xlim([solution.Time(range(1)) solution.Time(range(end)) ])
        ylim([min_v, max_v])
        ylabel('cpu time [s]')
        set(gca,'FontSize',15)
    hold off
end
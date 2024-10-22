
function f = make_plot_SOC(solution,Grid_para,range,modeColors)

set(groot,'defaultFigureVisible','on')

set(0, 'DefaultTextInterpreter', 'Latex')
set(0, 'DefaultLegendInterpreter', 'Latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'Latex')


f = figure('Renderer', 'painters', 'Position', [10 10 800 400]);
    hold on 
        data_soc= solution.SupCap_SOC;
        data_ref = 50*ones(length(solution.Time),1);
        
        min_v = -5 + min(min(data_soc,data_ref));
        max_v = 10 + max(max(data_soc,data_ref));
    
        plot_shade(solution,modeColors,min_v,max_v)
        plot_shade_name(solution,min_v,0.95*max_v)
    
        p1 = plot(solution.Time(range),data_soc(range),'LineWidth',2,'Color',[0, 0, 0]);
        p2 = plot(solution.Time(range),data_ref(range),'LineWidth',2,'LineStyle','--','Color',[0, 0, 0]);
        legend([p1 p2],'opf','state','Location','southwest','FontSize',15);
        xlim([solution.Time(range(1)) solution.Time(range(end)) ])
        ylim([min_v, max_v])
        ylabel('Super Cap SOC [%]')
        set(gca,'FontSize',15)
    hold off
end
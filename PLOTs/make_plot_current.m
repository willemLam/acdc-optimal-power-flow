
function f = make_plot_current(solution,Grid_para,constraints,range,modeColors)

set(groot,'defaultFigureVisible','on')

set(0, 'DefaultTextInterpreter', 'Latex')
set(0, 'DefaultLegendInterpreter', 'Latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'Latex')


f = figure('Renderer', 'painters', 'Position', [10 10 850 400]);
    hold on 
        data_opf= solution.reconstructed_Iabs(:,10)*Grid_para.I_b;
        data_state= abs(solution.variables_I(:,10))*Grid_para.I_b;
        data_ref = ceil(constraints.I_max(10)*ones(length(solution.Time),1)*Grid_para.I_b);
        
        min_v = -2 + min(min([data_opf(range),data_ref(range),data_state(range)]));
        max_v = 5 + max(max([data_opf(range),data_ref(range),data_state(range)]));
    
        plot_shade(solution,modeColors,min_v,max_v)
        plot_shade_name(solution,min_v,0.95*max_v)
    
        p1 = plot(solution.Time(range),data_opf(range),'LineWidth',2,'Color',[0.6, 0.6, 0.6]);
        p2 = plot(solution.Time(range),data_state(range),'LineWidth',2,'Color',[0, 0, 0]);
        p3 = plot(solution.Time(range),data_ref(range),'LineWidth',2,'LineStyle','--','Color',[0, 0, 0]);
        legend([p1 p2 p3],'opf','state','Ampacity','Location','southwest','FontSize',15);
        xlim([solution.Time(range(1)) solution.Time(range(end)) ])
        ylim([min_v, max_v])
        ylabel('Current B10-B11 [A]')
        set(gca,'FontSize',15)
    hold off
end
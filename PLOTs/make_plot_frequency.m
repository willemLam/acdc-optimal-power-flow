
function f = make_plot_frequency(solution,Grid_para,constraints,range,modeColors)

set(groot,'defaultFigureVisible','on')

set(0, 'DefaultTextInterpreter', 'Latex')
set(0, 'DefaultLegendInterpreter', 'Latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'Latex')

colorblindSet = {[215,25,28]/255;[253,174,97]/255;[171,217,233]*.8/255;[44,123,182]*.8/255};


f = figure('Renderer', 'painters', 'Position', [10 10 850 400]);
    hold on 
        data_Vup= abs(solution.f_upstream(:));
        data_Vdown = movmean(abs(solution.f_downstream(:)),5);
        
        min_v = -0.01 + min(min([data_Vup,data_Vdown]));
        max_v = 0.01 + max(max([data_Vup,data_Vdown]));
    
        plot_shade(solution,modeColors,min_v,max_v)
        plot_shade_name(solution,min_v,0.95*max_v)
    
        p1 = plot(solution.Time(range),data_Vup(range),'LineWidth',2,'Color',colorblindSet{2});
        p2 = plot(solution.Time(range),data_Vdown(range),'LineWidth',2,'Color',colorblindSet{4});
        legend([p1 p2],'Upstream','Downstream','Location','southwest','FontSize',15);
        xlim([solution.Time(range(1)) solution.Time(range(end)) ])
        ylim([min_v, max_v])
        ylabel('Frequency [Hz]','Interpreter','latex')
        set(gca,'FontSize',15)
    hold off
end
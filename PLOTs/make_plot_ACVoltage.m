
function f = make_plot_ACVoltage(solution,Grid_para,constraints,range,modeColors)

set(groot,'defaultFigureVisible','on')

set(0, 'DefaultTextInterpreter', 'Latex')
set(0, 'DefaultLegendInterpreter', 'Latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'Latex')

colorblindSet = {[215,25,28]/255;[253,174,97]/255;[171,217,233]*.8/255;[44,123,182]*.8/255};


f = figure('Renderer', 'painters', 'Position', [10 10 850 400]);
    hold on 
        data_V19= abs(solution.variables_E(:,19)*Grid_para.V_b);
        data_Vup= abs(solution.V_upstream*Grid_para.V_b);
        data_Vdown = abs(solution.V_downstream*Grid_para.V_b);
        
        min_v = -2 + min(min([data_V19(range),data_Vup(range),data_Vdown(range)]));
        max_v = 3 + max(max([data_V19(range),data_Vup(range),data_Vdown(range)]));
    
        plot_shade(solution,modeColors,min_v,max_v)
        plot_shade_name(solution,min_v,0.95*max_v)
    
        p1 = plot(solution.Time(range),data_V19(range),'LineWidth',2,'Color',[0.2 0.2 0.2]);
        p2 = plot(solution.Time(range),data_Vup(range),'LineWidth',2,'Color',colorblindSet{2});
        p3 = plot(solution.Time(range),data_Vdown(range),'LineWidth',2,'Color',colorblindSet{4});
        legend([p1 p2 p3],'Forming IC','Upstream','Downstream','Location','southwest','FontSize',15);
        xlim([solution.Time(range(1)) solution.Time(range(end)) ])
        ylim([min_v, max_v])
        ylabel('AC voltages [V]')
        set(gca,'FontSize',15)
    hold off
end
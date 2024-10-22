
function f = make_plot_DCVoltage(solution,Grid_para,constraints,range,modeColors)

set(groot,'defaultFigureVisible','on')

set(0, 'DefaultTextInterpreter', 'Latex')
set(0, 'DefaultLegendInterpreter', 'Latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'Latex')

colorblindSet = {[215,25,28]/255;[253,174,97]/255;[171,217,233]*.8/255;[44,123,182]*.8/255};


f = figure('Renderer', 'painters', 'Position', [10 10 850 400]);
    hold on 
        data_V22= abs(solution.variables_E(:,22)*Grid_para.Vdc_b);
        data_V23= abs(solution.variables_E(:,23)*Grid_para.Vdc_b);
        data_V24 = abs(solution.variables_E(:,24)*Grid_para.Vdc_b);
        
        min_v = -2 + min(min([data_V22(range),data_V23(range),data_V24(range)]));
        max_v = 3 + max(max([data_V22(range),data_V23(range),data_V24(range)]));
    
        plot_shade(solution,modeColors,min_v,max_v)
        plot_shade_name(solution,min_v,0.95*max_v)
    
        p1 = plot(solution.Time(range),data_V22(range),'LineWidth',2,'Color',[0.2 0.2 0.2]);
        p2 = plot(solution.Time(range),data_V23(range),'LineWidth',2,'Color',colorblindSet{2});
        p3 = plot(solution.Time(range),data_V24(range),'LineWidth',2,'Color',colorblindSet{4});
        legend([p1 p2 p3],'B22','B23','B24','Location','southwest','FontSize',15);
        xlim([solution.Time(range(1)) solution.Time(range(end)) ])
        ylim([min_v, max_v])
        ylabel('DC voltages [V]')
        set(gca,'FontSize',15)
    hold off
end
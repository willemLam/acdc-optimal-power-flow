
function f = make_plot_SOC_P(solution,Grid_para,range,modeColors)

set(groot,'defaultFigureVisible','on')

set(0, 'DefaultTextInterpreter', 'Latex')
set(0, 'DefaultLegendInterpreter', 'Latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'Latex')


f = figure('Renderer', 'painters', 'Position', [10 10 850 400]);
    hold on 
    %data for soc
        data_soc= solution.SupCap_SOC;
        data_ref = 50*ones(length(solution.Time),1);
    %data for P   
        data_P = solution.SupCap_P*Grid_para.Adc_b;

        delta_soc = max(abs(50-min([data_soc(range);data_ref(range)])),abs(50-max([data_soc(range);data_ref(range)])));
        min_soc = -5 + 50 - delta_soc;
        max_soc = 5 + 50 + delta_soc;

        delta_P = max(abs(min(data_P(range))),abs(max(data_P(range))));
        min_P = -500 - delta_P;
        max_P = 500 + delta_P;

        plot_shade(solution,modeColors,min_soc,max_soc)
        plot_shade_name(solution,min_soc,0.95*max_soc)
    
        % Left Y-Axis: Super Cap SOC
        yyaxis left
        p1 = plot(solution.Time(range),data_soc(range),'LineWidth',2,'Color',[0, 0, 0]);
        p2 = plot(solution.Time(range),data_ref(range),'LineWidth',2,'LineStyle','--','Color',[0, 0, 0]);
        
        ylabel('Super Cap SOC [\%]','Interpreter','latex')
        set(gca,'FontSize',15)
        ylim([min_soc, max_soc])

        % Right Y-Axis: Power
        yyaxis right
        p3 = plot(solution.Time(range), data_P(range), 'LineWidth', 2, 'Color', [0.6350, 0.0780, 0.1840]);  % Plot power in red
        ylabel('Power [W]','Interpreter','latex')  % Adjust units according to your data
        ylim([min_P, max_P])
        set(gca, 'YColor', [0.6350, 0.0780, 0.1840])

        % Set LaTeX interpreter explicitly for both y-axes
        ax = gca;
        ax.YAxis(1).Label.Interpreter = 'latex';
        ax.YAxis(2).Label.Interpreter = 'latex';

        legend([p1 p3], {'SOC', 'Power'}, 'Location', 'southwest', 'FontSize', 15)
        xlim([solution.Time(range(1)) solution.Time(range(end)) ])

        set(gca,'FontSize',15)
    hold off
end
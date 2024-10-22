
function f = make_plot_GCP(solution,Grid_para,constraints,range,modeColors)


set(groot,'defaultFigureVisible','on')

set(0, 'DefaultTextInterpreter', 'Latex')
set(0, 'DefaultLegendInterpreter', 'Latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'Latex')

colorblindSet = {[215,25,28]/255;[253,174,97]/255;[171,217,233]*.8/255;[44,123,182]*.8/255};



f = figure('Renderer', 'painters', 'Position', [10 10 850 400]);
    hold on 
    %data for soc
        data_P= real(solution.variables_S(:,1))*Grid_para.A_b;
    %data for P   
        data_Q = imag(solution.variables_S(:,1))*Grid_para.A_b;

        min_P = -500 + min(data_P(range));
        max_P = 500 + max(data_P(range));

        min_Q = -500 + min(data_Q(range));
        max_Q = 500 + max(data_Q(range));

        min_PQ = min(min_P,min_Q);
        max_PQ = max(max_P,max_Q);

        plot_shade(solution,modeColors,min_PQ,max_PQ)
        plot_shade_name(solution,min_PQ,0.95*max_PQ)
    
        % Left Y-Axis: Super Cap SOC
        yyaxis left
        p1 = plot(solution.Time(range),data_P(range),'LineWidth',2,'Color',[0, 0, 0]);

        ylabel('Active power [W]','Interpreter','latex')
        set(gca,'FontSize',15)
        ylim([min_PQ, max_PQ])

        % Right Y-Axis: Power
        yyaxis right
        p3 = plot(solution.Time(range), data_Q(range), 'LineWidth', 2, 'Color', colorblindSet{4});  % Plot power in red
        ylabel('Reactive power [VAr]','Interpreter','latex')  % Adjust units according to your data
        ylim([min_PQ, max_PQ])
        set(gca, 'YColor', colorblindSet{4})

        % Set LaTeX interpreter explicitly for both y-axes
        ax = gca;
        ax.YAxis(1).Label.Interpreter = 'latex';
        ax.YAxis(2).Label.Interpreter = 'latex';

        legend([p1 p3], {'Active', 'Reactive'}, 'Location', 'southwest', 'FontSize', 15)
        xlim([solution.Time(range(1)) solution.Time(range(end)) ])

        set(gca,'FontSize',15)
    hold off
end
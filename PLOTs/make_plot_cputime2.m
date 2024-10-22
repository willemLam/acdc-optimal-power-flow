
function f = make_plot_cputime2(solution,range)

%%
set(groot,'defaultFigureVisible','on')
colorgrd = 'blue_up';
m_cmap = colorgrad(8,colorgrd);
m_cmap2 = [0 0.4470 0.7410 ; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560;0.3010 0.7450 0.9330];


set(groot,'defaultFigureVisible','on')

set(0, 'DefaultTextInterpreter', 'Latex')
set(0, 'DefaultLegendInterpreter', 'Latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'Latex')

%%

f = figure('Renderer', 'painters', 'Position', [10 10 850 400]);

hold on
h1 = cdfplot(solution.duration(range));

set( h1, 'LineWidth', 2, 'Color', m_cmap(1,:));

% legend('Island OPF','Connected OPF','Location','SouthEast')
xlabel('Time [s]')
ylabel('Instances')
title('')

set(gca,'FontSize',20) 

end





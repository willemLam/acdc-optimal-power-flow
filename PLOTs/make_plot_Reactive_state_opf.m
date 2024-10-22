function f = make_plot_Reactive_state_opf(solution,Grid_para,range,modeColors)

f = figure('Renderer', 'painters', 'Position', [10 10 850 850]);
tiledlayout(4,1, 'TileSpacing', 'Compact'); 

m = 1;

nexttile
    hold on 
        data_opf = movmean(imag(solution.reconstructed_S(range,9))*Grid_para.A_b,m);
        data_state = movmean(imag(solution.variables_S(range,9))*Grid_para.A_b,m);
        
        min_v = -400+min(min(data_opf,data_state));
        max_v =  400+max(max(data_opf,data_state));
    
        plot_shade(solution,modeColors,min_v,max_v)
        plot_shade_name(solution,min_v,max_v)
    
        p2 = plot(solution.Time(range),data_state,'LineWidth',2,'Color',[0.6, 0.6, 0.6]);
        p1 = plot(solution.Time(range),data_opf,'LineWidth',2,'Color',[0, 0, 0]);
        
        legend([p1 p2],'opf','state','Location','southwest','FontSize',15);
        xlim([solution.Time(range(1)) solution.Time(range(end)) ])
        ylim([min_v, max_v])
        ylabel('BESS [VAr]')
        set(gca,'FontSize',15)
        set(gca,'xtick',[])
    hold off

nexttile
    hold on 
        data_opf = movmean(imag(solution.reconstructed_S(range,19))*Grid_para.A_b,m);
        data_state = movmean(imag(solution.variables_S(range,19))*Grid_para.A_b,m);
        
        min_v = -400+min(min(data_opf,data_state));
        max_v =  400+max(max(data_opf,data_state));
        plot_shade(solution,modeColors,min_v,max_v)
        
        p2 = plot(solution.Time(range),data_state,'LineWidth',2,'Color',[0.6, 0.6, 0.6]);
        p1 = plot(solution.Time(range),data_opf,'LineWidth',2,'Color',[0, 0, 0]);
    %     legend([p1 p2],'IC1-OPF','IC1-state','Location','southwest')
        xlim([solution.Time(range(1)) solution.Time(range(end)) ])
        ylim([min_v, max_v])
        ylabel('IC 1 [VAr]')
        set(gca,'FontSize',15)
        set(gca,'xtick',[])
    hold off

nexttile
    hold on 
        data_opf = movmean(imag(solution.reconstructed_S(range,20))*Grid_para.A_b,m);
        data_state = movmean(imag(solution.variables_S(range,20))*Grid_para.A_b,m);
        
        min_v = -400+min(min(data_opf,data_state));
        max_v =  400+max(max(data_opf,data_state));
        plot_shade(solution,modeColors,min_v,max_v)
        
        p2 = plot(solution.Time(range),data_state,'LineWidth',2,'Color',[0.6, 0.6, 0.6]);
        p1 = plot(solution.Time(range),data_opf,'LineWidth',2,'Color',[0, 0, 0]);
    %     legend([p1 p2],'IC2-OPF','IC2-state','Location','southwest')
        xlim([solution.Time(range(1)) solution.Time(range(end)) ])
        ylim([min_v, max_v])
        ylabel('IC 2  [VAr]')
        set(gca,'FontSize',15)
        set(gca,'xtick',[])
    hold off

 nexttile
    hold on 
        data_opf = movmean(imag(solution.reconstructed_S(range,21))*Grid_para.A_b,m);
        data_state = movmean(imag(solution.variables_S(range,21))*Grid_para.A_b,m);
        
        min_v = -400+min(min(data_opf,data_state));
        max_v =  400+max(max(data_opf,data_state));
        plot_shade(solution,modeColors,min_v,max_v)
        
        p2 = plot(solution.Time(range),data_state,'LineWidth',2,'Color',[0.6, 0.6, 0.6]);
        p1 = plot(solution.Time(range),data_opf,'LineWidth',2,'Color',[0, 0, 0]);
    %     legend([p1 p2],'IC3-OPF','IC3-state','Location','southwest')
        xlim([solution.Time(range(1)) solution.Time(range(end)) ])
        ylim([min_v, max_v])
        ylabel('IC 3 [VAr]')
        set(gca,'FontSize',15)
    hold off
end
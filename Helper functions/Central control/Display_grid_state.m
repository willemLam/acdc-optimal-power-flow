function Display_grid_state(solution,variable,Grid_para,idx)

    disp(['Current through line B10-B11: ', num2str(  round( solution.Iabs(10)*Grid_para.I_b)                              ,5), '  A']);
    disp(['Reactive power at slack: ', num2str(       round( solution.Q(1)*Grid_para.A_b)                          ,5), '  VAr']);
    disp(['DC voltage of ICs:       ', num2str(     round( solution.Eabs(idx.vscdc_vq)' *Grid_para.Vdc_b)            ,5), '  V']);
    disp(['Reactive power of ICs:     ', num2str(     round( solution.Q(idx.vscac_vq)' *Grid_para.A_b)          ,5), '  VAr']);
    disp(['Active power of ICs:       ', num2str(     round( solution.P(idx.vscac_vq)' *Grid_para.A_b)          ,5), '  W']);
   

end
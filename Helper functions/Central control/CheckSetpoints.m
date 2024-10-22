function proceed = CheckSetpoints(solution, constraints, Grid_para, idx)

    LF_threshold = 1e-3;
    Delta_E = max(abs(solution.Eabs - abs(solution.reconstructed.E)));
    Delta_Q = max(abs(solution.Q([idx.pqac; idx.vscac_vq]) - imag(solution.reconstructed.S([idx.pqac; idx.vscac_vq]))));
    Delta_P = max(abs(solution.P([idx.pqac; idx.pdc]) - real(solution.reconstructed.S([idx.pqac; idx.pdc]))));

 
    sp_check  = [(Delta_E<LF_threshold);
            (Delta_Q<LF_threshold);
            (Delta_P<LF_threshold);
            (solution.Eabs <= constraints.E_max);
            (solution.Eabs >= constraints.E_min);
            (solution.Iabs <= 1.1*constraints.I_max);
            (solution.Q(idx.vscac_vq) <= 0.2);
            (abs(solution.S(idx.vscac_vq)) <= 0.3)];

    if prod(sp_check)
        proceed = true;
        disp('The setpoints are valid')
    else
        fprintf(2,'Validate the setpoints \n')
        proceed = false;
    end

end
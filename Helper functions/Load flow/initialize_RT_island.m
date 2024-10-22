function [Grid_para, Filter_para, idx1, idx3, constraints] = initialize_RT_island()

    %% Base values
    A_b = 1e5;
    V_b= 400;
    Y_b = A_b/V_b^2; 
    I_b = A_b/(V_b*sqrt(3));
    
    Vdc_b = 800;
    Adc_b = A_b;
    Ydc_b = Adc_b/Vdc_b^2; 
    Idc_b = Adc_b/Vdc_b;
    
    %% Set the filter parameters
    Filter_para.IGBT_piecewise = [  0                   0
                                    0.04926559627563	0.7
                                    2.30625399327864	0.8
                                    15.7793399043317	0.85
                                    107.547461516782	0.8999
                                    735.837403888342	0.9499
                                    1588.01477341768	0.9699];
    Filter_para.Include_losses = 0;
    Filter_para.Exclude_losses = 1;
    Filter_para.R = 0; %checked
    Filter_para.X = 0;  %checked
    
    %% Set the Grid parameters
    Grid_para.n_dc = 8;
    Grid_para.n_ac = 22;
    Grid_para.n_ph = 1;
    Grid_para.n_nodes = Grid_para.n_ac*Grid_para.n_ph + Grid_para.n_dc;
    Grid_para.V_b = V_b;
    Grid_para.Y_b = Y_b;
    Grid_para.A_b = A_b;
    Grid_para.I_b = I_b;
    Grid_para.Vdc_b = Vdc_b;
    Grid_para.Ydc_b = Ydc_b;
    Grid_para.Adc_b = Adc_b;
    Grid_para.Idc_b = Idc_b;
    
    [Yac, YYLac, YLac, YTac, YYTac, I_bac, Ampacitiesac, y_lx, y_tx, Aac, linedata_ac]  = Ymatrix_non_sym_shunt('linedata_AC_LCL_island2.txt',A_b,V_b,[]);
    [Ydc, YYLdc, YLdc, YTdc, YYTdc, I_bdc, Ampacitiesdc, y_ihdc, y_idc, Adc, linedata_dc]  = Ymatrix_non_sym_shunt('linedata_DC_LCL.txt',Adc_b,Vdc_b,[]);
    Ydc = Ydc(end-Grid_para.n_dc + 1:end,end-Grid_para.n_dc + 1:end)*2;
    Yac = cell2mat(arrayfun(@(x) x*eye(Grid_para.n_ph),Yac,'UniformOutput',false));
    YY = blkdiag(Yac,Ydc);
    
    YYLdc = YYLdc(end-Grid_para.n_dc + 1:end,end-Grid_para.n_dc + 1:end)*2;
    YYLac = cell2mat(arrayfun(@(x) x*eye(Grid_para.n_ph),YYLac,'UniformOutput',false));
    YYL = blkdiag(YYLac,YYLdc);
    
    YYTdc = YYTdc(end-Grid_para.n_dc + 1:end,end-Grid_para.n_dc + 1:end)*2;
    YYTac = cell2mat(arrayfun(@(x) x*eye(Grid_para.n_ph),YYTac,'UniformOutput',false));
    YYT = blkdiag(YYTac,YYTdc);
    
    Ampacitiesdc = Ampacitiesdc(end-Grid_para.n_dc + 1:end,end-Grid_para.n_dc + 1:end);
    Ampacitiesac = cell2mat(arrayfun(@(x) x*eye(Grid_para.n_ph),Ampacitiesac,'UniformOutput',false));
    Ampacities = blkdiag(Ampacitiesac,Ampacitiesdc);
    
    A = blkdiag(Aac,Adc(:,end-7:end));
    Grid_para.E_0 = [ones(Grid_para.n_ac,1) ; ones(Grid_para.n_dc,1)];
    
    Grid_para.G = real(YY);
    Grid_para.B = imag(YY);
    Grid_para.Ampacities = Ampacities;
    Grid_para.YY = YY;
    Grid_para.YYL = YYL;
    Grid_para.YYT = YYT;
    Grid_para.Yac = Yac;
    Grid_para.Ydc = Ydc;
    Grid_para.A = A;
    Grid_para.n_lines = size(linedata_ac,1) + size(linedata_dc,1);

    Grid_para.lines = [linedata_ac(:,1:2);linedata_dc(:,1:2)];
    Grid_para.I_indexu = sub2ind(size(Ampacities),Grid_para.lines(:,1),Grid_para.lines(:,2));
    Grid_para.Ampacities_list = Ampacities(Grid_para.I_indexu);
    
    %% Set the nodes types
    idx1.slack = [];
    idx1.pqac = [1:18]';
    idx1.pvac = []';
    
    idx1.pdc = [28:30]';
    idx1.vdc = [27]';
    
    idx1.vscac_pq = []';
    idx1.vscac_vq = [20:22]';
    idx1.vscac_vv = [19]';
    
    idx1.vscdc_pq = []';
    idx1.vscdc_vq = [24:26]';
    idx1.vscdc_vv = [23]';
    
    idx3 = Get_multiphase_Node_indices(idx1,Grid_para);
    linedata = [linedata_ac;linedata_dc];
    Grid_para = Get_Converter_para(idx1,linedata,Grid_para);

%     disp('Initialization RT control');

    %% Data from SCADA
    Grid_para.SCADA_bus = [1:30];
    Grid_para.SCADA_line = [1:21,22:26,26,26];


    %% Constraints
    constraints.IC1_Qmin = -0.3;
    constraints.IC1_Qmax =  0.3;
    constraints.IC2_Qmin = -0.3;
    constraints.IC2_Qmax =  0.3;
    constraints.IC3_Qmin = -0.3;
    constraints.IC3_Qmax =  0.3;
    constraints.IC4_Qmin = -0.3;
    constraints.IC4_Qmax =  0.3;

    constraints.PV_Q_perun_min = -0.3;
    constraints.PV_Q_perun_max =  0.3;

    constraints.E_min = [0.95*ones(Grid_para.n_ac,1); 0.900*ones(Grid_para.n_dc,1)];
    constraints.E_max = [1.05*ones(Grid_para.n_ac,1); 0.975*ones(Grid_para.n_dc,1)];

    constraints.I_max = Ampacities(Grid_para.I_indexu);

    constraints.Sambat_Pmax = 0.3;
    constraints.Sambat_Pmin = -0.3;

    constraints.SupCap_Pmax =  0.05;
    constraints.SupCap_Pmin = -0.05;

    if Grid_para.n_ph == 3
        Grid_para.E_0 = [repmat([1; exp(-1i*pi*2/3);  exp(1i*pi*2/3) ], Grid_para.n_ac,1 ) ; ones(Grid_para.n_dc,1)];
    elseif Grid_para.n_ph == 1
        Grid_para.E_0 = [repmat([1], Grid_para.n_ac,1 ) ; ones(Grid_para.n_dc,1)];
    end

end
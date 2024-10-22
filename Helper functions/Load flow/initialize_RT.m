function [Grid_para, Filter_para, idx1, idx3, constraints] = initialize_RT(topology)

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
    
    Filter_para.Include_losses = 1;
    % b* Iac + c* Iac^2 + a +  d* Idc + e* Idc^2
%     Filter_para.a = [ 0.0045     0.0045    0.0045  ]'; %pu
%     Filter_para.b = [ 0         0        0      ]';
%     Filter_para.c = [ 0         0        0      ]';
%     Filter_para.d = [ 0.01       0.02    0.02    ]';
%     Filter_para.e = [ 0.1         0.2      0.2      ]';
    
%     Filter_para.b = [ 0         0        0     ]';
%     Filter_para.c = [ 0.0028    0.0028   0.0028]';
%     Filter_para.a = [ 0.0053    0.0053   0.0053]'; %pu
%     Filter_para.d = [ 0.1120    0.1120   0.1120]';
%     Filter_para.e = [-0.9632   -0.9632  -0.9632]';

    Filter_para.b = [ 0         0        0     ]';
    Filter_para.c = [ 0.0028    0.0028   0.0028]';
    Filter_para.a = [ 0.0059    0.0059   0.0059]'; %pu
    Filter_para.d = [ 0.0903    0.0903   0.0903]';
    Filter_para.e = [ 0         0        0]';

%     Filter_para.a = -[ 0.0001    -0.008          -0.008]'; %pu
%     Filter_para.b = -[ 0         0              0]';
%     Filter_para.c = -[-0.5498    0.3020        0.3020]';
%     Filter_para.d = -[-0.1076    0.0884        0.0884]';
%     Filter_para.e = -[ 1.5865    -1.2769         -1.2769]';

%     % a + b* Iac + c* Iac^2 + d* Idc + e* Idc^2
%     Filter_para.a = [ 0.0025   -0.0107         0]'; %pu
%     Filter_para.b = [-0.0685    0.1053         0]';
%     Filter_para.c = [-0.1302   -0.3940         0]';
%     Filter_para.d = [-0.1020    0.0627         0]';
%     Filter_para.e = [ 1.4006   -1.0241         0]';

    %% Set the Grid parameters
    Grid_para.n_dc = 6;
    Grid_para.n_ac = 21;
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
    
    [Yac, YYLac, YLac, YTac, YYTac, I_bac, Ampacitiesac, y_lx, y_tx, Aac, linedata_ac]  = Ymatrix_non_sym_shunt('linedata_AC_LCL_island.txt',A_b,V_b,[]);
    [Ydc, YYLdc, YLdc, YTdc, YYTdc, I_bdc, Ampacitiesdc, y_ihdc, y_idc, Adc, linedata_dc]  = Ymatrix_non_sym_shunt('linedata_DC_LCL_island.txt',Adc_b,Vdc_b,[]);
    Ydc = Ydc(end-Grid_para.n_dc + 1:end,end-Grid_para.n_dc + 1:end)*1; %*2
    Yac = cell2mat(arrayfun(@(x) x*eye(Grid_para.n_ph),Yac,'UniformOutput',false));
    YY = blkdiag(Yac,Ydc);
    
    YYLdc = YYLdc(end-Grid_para.n_dc + 1:end,end-Grid_para.n_dc + 1:end)*1; %*2
    YYLac = cell2mat(arrayfun(@(x) x*eye(Grid_para.n_ph),YYLac,'UniformOutput',false));
    YYL = blkdiag(YYLac,YYLdc);
    
    YYTdc = YYTdc(end-Grid_para.n_dc + 1:end,end-Grid_para.n_dc + 1:end)*1; %*2
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
    if topology == "CONNECTED"
        idx1.slack = 1;
        idx1.pqac = [2:18]';
        idx1.pvac = []';
        
        idx1.pdc = [25:27]';
        idx1.vdc = []';     %27
        
        idx1.vscac_pq = [20:21]';
        idx1.vscac_vq = [19]';
        idx1.vscac_vv = []';
        
        idx1.vscdc_pq = [23:24]';
        idx1.vscdc_vq = [22]';
        idx1.vscdc_vv = []';
    
        idx3 = Get_multiphase_Node_indices(idx1,Grid_para);
        linedata = [linedata_ac;linedata_dc];
        Grid_para = Get_Converter_para(idx1,linedata,Grid_para);

    elseif topology == "ISLAND"
        idx1.slack = [];
        idx1.pqac = [1:18]';
        idx1.pvac = []';
        
        idx1.pdc = [26:27]';
        idx1.vdc = [25]';
        
        idx1.vscac_pq = [20:21]';
        idx1.vscac_vq = []';
        idx1.vscac_vv = [19]';
        
        idx1.vscdc_pq = [23:24]';
        idx1.vscdc_vq = []';
        idx1.vscdc_vv = [22]';
    
        idx3 = Get_multiphase_Node_indices(idx1,Grid_para);
        linedata = [linedata_ac;linedata_dc];
        Grid_para = Get_Converter_para(idx1,linedata,Grid_para);

    else
        error()
    end
   

    %% Data from SCADA
    Grid_para.SCADA_bus = [1:21,23:25,27:29];
    Grid_para.SCADA_line = [1:20,22:24,26,27];


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

    constraints.E_min = [0.95*ones(Grid_para.n_ac,1); 0.875*ones(Grid_para.n_dc,1)];
    constraints.E_max = [1.05*ones(Grid_para.n_ac,1); 0.975*ones(Grid_para.n_dc,1)];

    constraints.I_max = Grid_para.Ampacities_list;

    constraints.Sambat_Pmax =  0.25; %0.3;
    constraints.Sambat_Pmin = -0.25; %-0.3;

    constraints.SupCap_Pmax =  0.03; %0.05;
    constraints.SupCap_Pmin = -0.03; %-0.05;

    constraints.BESS_loss = 355/Grid_para.A_b; 

    if Grid_para.n_ph == 3
        Grid_para.E_0 = [repmat([1; exp(-1i*pi*2/3);  exp(1i*pi*2/3) ], Grid_para.n_ac,1 ) ; ones(Grid_para.n_dc,1)];
    elseif Grid_para.n_ph == 1
        Grid_para.E_0 = [repmat([1], Grid_para.n_ac,1 ) ; ones(Grid_para.n_dc,1)];
    end


    constraints.Edc_nominal = 720/800; % variables.E_slack_sp; %

end
function [E_bus, S_bus, I_line] = Get_SE_data(TYPE,UDP_port, Grid_para,i)

count = 0;
try
    
    if TYPE == 'SIMULATION'
%         S_bus = S_star;
%         S_bus(11) = P_max_solarmax+P_max_facade;
%         S_bus(14) = complex(-0.2,15000*sin((i)/100)/A_b);
%         if i == 1
%         S_bus(1) = -sum(S_bus);
%         end
%         [E_bus,~,~] = NR_rectangularACDC_1ph_general(Grid_para,Filter_para,S_bus,E,E_0,idx1,1e-7,20);
%         I_long = abs(get_Current_flow(E,Grid_para));
%         I_line = I_long(I_indexu);

%get saved state
        FileName_states = 'C:\Users\admin\Documents\git-clones\States\State_09-12-2023 17-00.mat';
        FileName_states = 'results\09-14-2023 16-32_EXPERIMENT_YES.mat';
        nameD = ['data_',num2str(i)];
        D2 = load(FileName_states, genvarname(nameD));
        d2= D2.(genvarname(nameD));
        
%get state in absolute units
%         E_bus = [d.E_scada(1:18)*(V_b/sqrt(3)); d.E_scada(19:26)*Vdc_b];
%         I_line = [d.I_scada(1:end-4)*(I_b); d.I_scada(end-3:end)*Idc_b ];
%         S_bus = d.S_scada/3*A_b;
% %get state in pu
%         S_bus = [S_bus(1:Grid_para.n_ac*Grid_para.n_ph,:)*3; S_bus(Grid_para.n_ac*Grid_para.n_ph+1:end,:)]/A_b;
% 
%         Eac_bus = E_bus(1:Grid_para.n_ac*Grid_para.n_ph,:)/(V_b/sqrt(3));
%         Edc_bus = E_bus(Grid_para.n_ac*Grid_para.n_ph+1:end,:)/Vdc_b;
%         E_bus = [Eac_bus; Edc_bus];
% 
%         I_long = abs(get_Current_flow(E_bus,Grid_para));
%         I_line = I_long(I_indexu);

        E_bus = d2.E_scada;
        I_line = d2.I_scada;
        S_bus = d2.S_scada;

        if length(E_bus) < Grid_para.n_nodes
            E_bus = [E_bus(1:22); E_bus(19:end) ];
            S_bus = [S_bus(1:14); zeros(4,1); S_bus(15:end) ];
            I_line = [I_line(1:17); I_line(14:end) ];
        end

    elseif TYPE == 'EXPERIMENT'

        E_all = []; S_all = [];  I_all = [];

        for j = 1:5
        [E_all(:,j),S_all(:,j),I_all(:,j)] = get_states_from_SE(UDP_port, Grid_para);
        end

        E_bus = median(E_all,2);
        S_bus = median(S_all,2);
        I_line = median(I_all,2);
        S_bus(19:end) = -S_bus(19:end);
    end

   

    E_bus(abs(E_bus)<0.5) = 0.9375;%   %numerical trick to keep computing the SC when the DC grid is off

%     [E_bus_rec,~,~] = NR_rectangularACDC_1ph_general(Grid_para,Filter_para,S_bus,E_bus,E_0,idx1,1e-7,20) 
%     E_bus_rec.*conj(YY*E_bus_rec)


catch
    disp('Cannot connect to scada')
    count = count + 1;
    if count == 3
        error('Check scada connection')
    end
end    

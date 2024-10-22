
function [E_bus,S_bus,I_line] = get_states_from_SE(UDP_port, Grid_para)

n_ph = Grid_para.n_ph;
n_dc = Grid_para.n_dc;
n_ac = Grid_para.n_ac;
n_nodes = Grid_para.n_nodes;
V_b = Grid_para.V_b;
Y_b = Grid_para.Y_b;
A_b = Grid_para.A_b;
Vdc_b = Grid_para.Vdc_b;
Ydc_b = Grid_para.Ydc_b;
Adc_b = Grid_para.Adc_b;

I_b = A_b/(V_b*sqrt(3));
Idc_b = Adc_b/Vdc_b;

%% Get UDP


    [message, ~] = judp('RECEIVE',35005,2^18,4000); 
    json1 = char(message);
%     fprintf(1,'%s\n',json1);
    data = loadjson(json1');



    numBuses = size(data.Buses,2);
    numLines = size(data.Lines,2);
%% Extract P, Q, V, I
    for i = 1:numBuses
        for ph = 1:n_ph
            P_bus(i,ph) = data.Buses{1,i}.Terminals{1,ph}.P;  
            Q_bus(i,ph) = data.Buses{1,i}.Terminals{1,ph}.Q;
            V_bus_abs(i,ph) = data.Buses{1,i}.Terminals{1,ph}.V_abs;
            V_bus_arg(i,ph) = data.Buses{1,i}.Terminals{1,ph}.V_arg;
        end
    end



    for j = 1:numLines
        for ph = 1:n_ph
            I_line_abs(j,ph) = data.Lines{1,j}.Conductors{1,ph}.I_abs;
            I_line_arg(j,ph) = data.Lines{1,j}.Conductors{1,ph}.I_arg;
        end
    end



%% Transform to base
P_bus = [P_bus(Grid_para.SCADA_bus(1:n_ac))*3 ; P_bus(Grid_para.SCADA_bus(n_ac+4+1:end))];
Q_bus = [Q_bus(Grid_para.SCADA_bus(1:n_ac))*3 ; Q_bus(Grid_para.SCADA_bus(n_ac+4+1:end))];
S_bus = complex(P_bus,Q_bus)/A_b;

Eac_bus = V_bus_abs(Grid_para.SCADA_bus(1:n_ac))/(V_b/sqrt(3)) .* exp(1i*V_bus_arg(Grid_para.SCADA_bus(1:n_ac)));

% % little trick to avoid a lot of math
% Eangle_slack = imag(Eac_bus(19));
% Eac_bus = abs(Eac_bus).*(cos(angle(Eac_bus) - Eangle_slack) + 1i*sin(angle(Eac_bus) - Eangle_slack));

Edc_bus = V_bus_abs(Grid_para.SCADA_bus(n_ac+4+1:end))/Vdc_b;
E_bus = [Eac_bus; Edc_bus];

Iac_line = I_line_abs(Grid_para.SCADA_line(1:21))/(I_b) .* exp(1i*I_line_arg(Grid_para.SCADA_line(1:21)));
Idc_line = I_line_abs(Grid_para.SCADA_line(22:end))/Idc_b;
I_line = [Iac_line; Idc_line];

%% Select the LV and DC states
% 
% E_bus = E_bus(Grid_para.SCADA_bus);
% S_bus = S_bus(Grid_para.SCADA_bus);
% I_line = I_line(Grid_para.SCADA_line);

end
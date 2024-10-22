function variables = Receive_data_from_SCADA(Grid_para,variables,constraints,idx,topology)

% info


% CDA_port = 35211;
% CDA_host = '192.168.1.86';

SCADA_port = 35005;
SCADA_host = '192.168.1.2';

count = 0;

try

[message, ~] = judp('RECEIVE',SCADA_port,2^15,4000); 
json1 = char(message);
data.SCADA = loadjson(json1');

% data = jsondecode(json1');



%  initialisation
alpha = exp(1i*pi*2/3);

T = 1/3* [1    1       1;
          1    alpha   alpha^2;
          1    alpha^2 alpha];

E = zeros(Grid_para.n_nodes,1);
S = zeros(Grid_para.n_nodes,1);
I = zeros(Grid_para.n_lines,1);

            %% ----------- %%
            %%    SCADA    %%
            %% ----------- %%

%% AC network
for i = 1:22
    %voltage
    E_3ph = [data.SCADA.Data.Buses{1, i}.Terminals{1, 1}.V_abs  * exp(1i*data.SCADA.Data.Buses{1, i}.Terminals{1, 1}.V_arg);
             data.SCADA.Data.Buses{1, i}.Terminals{1, 2}.V_abs  * exp(1i*data.SCADA.Data.Buses{1, i}.Terminals{1, 2}.V_arg);
             data.SCADA.Data.Buses{1, i}.Terminals{1, 3}.V_abs  * exp(1i*data.SCADA.Data.Buses{1, i}.Terminals{1, 3}.V_arg)];
    E(i) = T(2,:)*E_3ph /(Grid_para.V_b/sqrt(3));

    %power
    S_3ph = [complex(data.SCADA.Data.Buses{1, i}.Terminals{1, 1}.P, data.SCADA.Data.Buses{1, i}.Terminals{1, 1}.Q);
             complex(data.SCADA.Data.Buses{1, i}.Terminals{1, 2}.P, data.SCADA.Data.Buses{1, i}.Terminals{1, 2}.Q);
             complex(data.SCADA.Data.Buses{1, i}.Terminals{1, 3}.P, data.SCADA.Data.Buses{1, i}.Terminals{1, 3}.Q)];
    S(i) = sum(S_3ph) / Grid_para.A_b;

end

%% DC network
for i = 23:30
    %voltage
    E(i) = data.SCADA.Data.Buses{1, i}.Terminals{1, 1}.V_abs /Grid_para.Vdc_b;
        if E(i) < 0.5 %little trick to text when DC is off
            E(i) = constraints.Edc_nominal;
        end
    %power
    S(i) = data.SCADA.Data.Buses{1, i}.Terminals{1, 1}.P /Grid_para.A_b;

end

%% AC network
for i = 1:21
    %current
    I_3ph = [data.SCADA.Data.Lines{1, i}.Conductors{1, 1}.I_abs  * exp(1i*data.SCADA.Data.Lines{1, i}.Conductors{1, 1}.I_arg);
             data.SCADA.Data.Lines{1, i}.Conductors{1, 2}.I_abs  * exp(1i*data.SCADA.Data.Lines{1, i}.Conductors{1, 2}.I_arg);
             data.SCADA.Data.Lines{1, i}.Conductors{1, 3}.I_abs  * exp(1i*data.SCADA.Data.Lines{1, i}.Conductors{1, 3}.I_arg)];
    I(i) = T(2,:)*I_3ph /(Grid_para.I_b);

end

%% DC network
for i = 22:27
    %current
    I(i) = data.SCADA.Data.Lines{1, i}.Conductors{1, 1}.I_abs /Grid_para.Idc_b;
end


E([22,26,30]) = [];
S([22,26,30]) = [];
I([21,25]) = [];

S = E.*conj(Grid_para.YY*E);
I = get_Current_flow(E,Grid_para); %we need the flow

%% add artificial PV to test constraint violation
% 
% S(11) =  0.15  ;
% 
% tol = 1e-7; n_max = 100;
% [~, Filter_para, ~, ~, ~] = initialize_RT("ISLAND");
% 
% [E,~,~,~] = NR_rectangularACDC_1ph_general_V2_quadratic_loss(Grid_para,Filter_para,S,E,E,idx,tol,n_max);
% S = E.*conj(Grid_para.YY*E);
% I = get_Current_flow(E,Grid_para);
%% 

variables.E = E;
variables.Eabs = abs(E);
variables.Eang = angle(E);

if topology == "ISLAND" %change the ref angle of the ac nodes
    variables.Eang(1:Grid_para.n_ac) = variables.Eang(1:Grid_para.n_ac) - variables.Eang(idx.vscac_vv);
end

variables.E = variables.Eabs.*exp(1i*variables.Eang);

% variables.S =  variables.E .* conj(Grid_para.YY*variables.E);
% variables.I = get_Current_flow(variables.E,Grid_para);


variables.S = S;
variables.P = real(S);
variables.Q = imag(S);

% 
variables.I = I;
variables.Iabs = abs(I);



%% Point of Common Coupling
V_downstream = T(2,:) * transpose((data.SCADA.Data.downstream.V_abs .* exp(1i * data.SCADA.Data.downstream.V_arg)));
V_upstream = T(2,:) * transpose((data.SCADA.Data.upstream.V_abs .* exp(1i * data.SCADA.Data.upstream.V_arg)));
f_downstream = data.SCADA.Data.downstream.frequency;
f_upstream = data.SCADA.Data.upstream.frequency;

variables.V_downstream = V_downstream/(Grid_para.V_b/sqrt(3));
variables.V_upstream   = V_upstream/(Grid_para.V_b/sqrt(3));
variables.f_downstream = f_downstream;
variables.f_upstream   = f_upstream;



catch
    count = count + 1;
    if count == 3
           warning('Cannot connect to SCADA')
    end
end


end
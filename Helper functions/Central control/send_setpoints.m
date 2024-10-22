function [] = send_setpoints(solution,variable,Grid_para,idx)

% solarmax
%         count = 0;
%         dest_port = 12345;
%         dest_IP = '192.168.1.132';
%         P = solution.P(11)*Grid_para.A_b/1e3 ;
%         Q = 0*Grid_para.A_b/1e3 ;
%         try
%         if P > 1 %converter will switch off
%            sendToSolarMax(P, Q, dest_port, dest_IP)
%         end
%         catch
%             count = count + 1;
%             if count == 3
%                    warning('Cannot connect to Interfacing converters')
%             end
%         end

        % Interfacing Converters
        count = 0;
        dest_port = 34691;
        dest_IP = "192.168.1.20";
        sp_IC_Q = solution.Q(idx.vscac_vq)*Grid_para.A_b; 
        sp_IC_Q(abs(sp_IC_Q) < 100) = 0; % no need to inject values <100
        sp_IC_E = solution.Eabs(idx.vscdc_vq)*Grid_para.Vdc_b;

        try
           % send in W
        sendToInterfacingConverter(sp_IC_Q, sp_IC_E, time, iter, dest_port, dest_IP)
        catch
            count = count + 1;
            if count == 3
                   warning('Cannot connect to Interfacing converters')
            end
        end
end
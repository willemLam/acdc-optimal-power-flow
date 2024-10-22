function [P_resource_IC,Q_resource_IC,Edc_resource_IC] = Get_IC_data(Grid_para)

    A_b = Grid_para.A_b;
    Vdc_b = Grid_para.Vdc_b;
    
    
    try
    [time_IC,P_resource_IC,Q_resource_IC,Edc_resource_IC] = readFromInterfacingConverter();
    P_resource_IC = P_resource_IC/A_b;
    Q_resource_IC = Q_resource_IC/A_b;
    Edc_resource_IC = Edc_resource_IC/Vdc_b;
    
    catch
        P_resource_IC = 0; %previous PV (accounting for the curtailment)
        Q_resource_IC = 0; %previous PV 
        Edc_resource_IC = 0; %previous PV 
        disp('Cannot read from IC - Previous IC data is used')

    end

end
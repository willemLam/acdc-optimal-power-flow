function [] =  sendToSolarMax(P, Q, dest_port, dest_IP)
% write to a specif host and port, with the
% first 4 arguments the values to be send in json format.
% P, Q in (kW, kVar)
% defined ports and IP for the resources 

% port = 12345;
% host = '192.168.1.132';
% P = 16; 
% 
port = dest_port;
host = dest_IP;
% -------------------------------------------------------------------------
%connection settings
number_of_retries = 1; % set to -1 for infinite
samplingrate=1000; %in ms
i = 0; 

%% set-points as percentage of rating
ratings = [8190 7980]; %[T1 T2] Watts

P_T1_perc = 1e3*P*ratings(1)/sum(ratings)*100/10e3;
P_T2_perc = 1e3*P*ratings(2)/sum(ratings)*100/10e3;

%% string for sending to T1

string_T1 = 'FB;69;1A|C8:PLR=Setpoint|';
P_T1_hex = dec2hex(floor(P_T1_perc));

if numel(P_T1_hex)<3
    for i =1:3- numel(P_T1_hex)
    P_T1_hex = strcat('0',P_T1_hex );
    end
end


string_T1 = strrep(string_T1,'Setpoint' , P_T1_hex);
append_T1 = dec2hex(sum(uint16(unicode2native(string_T1))));

if numel(append_T1)<4
    for i =1:4- numel(append_T1)
    append_T1 = strcat('0',append_T1 );
    end
end
    

string_T1 = strcat('{', string_T1, append_T1, '}');

%% string for sending to T2

string_T2 = 'FB;37;1A|C8:PLR=Setpoint|';
P_T2_hex = dec2hex(floor(P_T2_perc));

if numel(P_T2_hex)<3
    for j =1:3-numel(P_T2_hex)
    P_T2_hex = strcat('0',P_T2_hex );
    end
end

string_T2 = strrep(string_T2,'Setpoint' , P_T2_hex);
append_T2 = dec2hex(sum(uint16(unicode2native(string_T2))));

if numel(append_T2)<4
    for j =1:4-numel(append_T2)
    append_T2 = strcat('0',append_T2 );
    end
end
string_T2 = strcat('{', string_T2, append_T2, '}');

fprintf(['** sending setpoint solarmax **\n', 'P: ' num2str(P*1e3), ' Q: ' num2str(Q)], '\n');

%% sending
i = 0;
while(i == 0)
    disp('sending')
%     data_to_send=[regexprep(num2str(15000),'','')]
%     data_to_send = '{FB;69;1A|C8:PLR=002|0549}' %200
    % {FB;69;1A|C8:PLR=032|054C} - 105, 5000
    % {FB;37;1A|C8:PLR=032|0547} - 55, 5000
    % {FB;37;1A|C8:PLR=050|0547} - 55, 8000
%     data_to_send = '{FB;69;1A|C8:PLR=014|054C}'; % 2000
%     data_to_send = '{FB;69;2E|64:PAC;PDC;UD01;ID01;ID02;UD02|0A24}';
    solarMax_T1=tcpSendRecieve(port,host,number_of_retries,samplingrate,string_T1);
    solarMax_T2=tcpSendRecieve(port,host,number_of_retries,samplingrate,string_T2);

%     if(numel(data_received)>0 )
%         % Do your stuff with the data
% %         mean_data = mean(data_received)
% %         disp('the average value of data recieved=')
% %     disp(mean_data)
%     end;
    i = i+1;      
end


% 
% 
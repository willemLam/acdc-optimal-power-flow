% http://128.179.34.35:52100/d/EPUyyHw7k/microgrid-mpc-experiment?orgId=1

clear all

%% Get data

GET_DATA = false;
if GET_DATA
database = 'logging';

seriesname = 'Microgrid';
tag_keys = {'MPC'};
tag_values = {'P_load'};

t0 = convertTo(datetime(2019,9,10,0,0,0),'epochtime') - 3600*2
tf = convertTo(datetime(2019,9,10,23,59,59),'epochtime') - 3600*2

[P_load, ~] = readSeriesFromDatabase_InfluxDB2(database, seriesname, tag_keys, tag_values, [], t0, tf, 600);

Q_load = 0.7*P_load + randn(length(P_load),1)
figure
hold on
plot(interp(P_load,60))
plot(interp(Q_load,60))

S_load = complex(interp(P_load,60),interp(Q_load,60))
save('Zenone_profile.mat',"S_load")

end

Ts = 2;
load('Zenone_profile.mat',"S_load"); %in kW, kVAr
S_load = S_load/2;
t0_load = datetime('now')
[ h , m , s ] = hms(t0_load);
t_load = floor((s + m*60 + h*3600)/10);


for i = 1:2888
    time = now();
    iter = i;

    CC_full_start = tic;

        count = 0;
        dest_port = 34009;
        dest_IP = '192.168.1.18';
        P = -real(S_load(t_load + i))*1e3;
        Q = -imag(S_load(t_load + i))*1e3;
        try
           sendToZenone(P, Q, time, iter, dest_port, dest_IP)
        catch
            count = count + 1;
            if count == 3
                   warning('Cannot connect to Zenone')
            end
        end


    CC_full_time = toc(CC_full_start);
    time_to_wait = Ts-CC_full_time;
    pause(time_to_wait);


end
clear all
% library for JSON
addpath('jsonlab-1.5\')
addpath('judp\judp.m')

addpath(genpath(pwd))

%t = datenum(2016, 11, 18, 1, 1, 1);
t = now() - 0;
TIMEZONE = 1;
cur_second_of_day = round((t - today())*24*3600) - TIMEZONE*3600;
%cur_minute_of_day = round((t - tmp_day)*24*3600) - TIMEZONE*3600;
k_cur = floor(cur_second_of_day/60);
% k_cur = 800;
n_ph = 1;
for i = 1:10

    
UDP_port = 35000;
[E_bus,S_bus,I_line] = get_states_from_SE(UDP_port, Grid_para);


x2 = [1 2];
y2 = [4 4];
p = histogram(x2);
xlim([0 100])
ylim([2.5 4])
xlabel('Iteration')
ylabel('Approximation for \pi')

p.XDataSource = 'x2';
p.YDataSource = 'y2';

denom = 1;
k = -1;
for t = 3:100
    denom = denom + 2;
    x2(t) = t;
    y2(t) = 4*(y2(t-1)/4 + k/denom);
    refreshdata
    drawnow
    k = -k;
end

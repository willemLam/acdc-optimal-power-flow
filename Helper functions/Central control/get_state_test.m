FileName = fullfile(['C:\Users\admin\Documents\Willem\acdc_centralized_controller\test\',datestr(now,'mm-dd-yyyy HH-MM') ,'.mat']);
emptyStruct = struct;
save(FileName,'-struct','emptyStruct');
% m = matfile(FileName, 'Writable', true); %Note: writable is true by default IF the file does not exist

[Grid_para, Filter_para, idx1, idx3, constraints] = initialize_RT();
idx = idx1;
port_SE = 35005;


for i = 1:10000

    i;
    [message, ~] = judp('RECEIVE',port_SE,2^18,4000); 
    json1 = char(message);
%     fprintf(1,'%s\n',json1);
    data = loadjson(json1');

    solution.z = data.z;

    try
    solution.z(241)*800

    end
    name = ['data_',num2str(i)];
    assignin('base',name, solution)
    save(FileName,name,'-append')
    clear(genvarname(name))   

        
    
end

%%

clear z

% FileName = 'C:\Users\admin\Documents\Willem\acdc_centralized_controller\test\03-15-2024 14-54.mat';
% FileName = 'C:\Users\admin\Documents\Willem\acdc_centralized_controller\test\03-15-2024 15-46.mat';

%forming
%FileName = 'C:\Users\admin\Documents\Willem\acdc_centralized_controller\test\03-20-2024 11-01.mat'
%FileName = 'C:\Users\admin\Documents\Willem\acdc_centralized_controller\test\03-20-2024 11-15.mat'

s = whos('-file',FileName);
N = length({s.name});

z = [];

for j = 2:N
    j
name = ['data_',num2str(j)];
D = load(FileName, genvarname(name));
d= D.(genvarname(name));

    try
        z(:,j) = d.z;
    catch
        z(:,j) = z(:,j-1);
    end
end
   
%%

figure
subplot(4,1,1)
plot(800*real(z(241:248,:))')
ylabel('DC Voltage')


subplot(4,1,2)
plot(125*real(z(249:end,:))')
ylabel('DC Currents')

subplot(4,1,3)
plot(abs(complex(z_new(190:192,:),z_new(223:225,:)))')
ylabel('AC Currents')

subplot(4,1,4)
plot(abs(complex(z(187:189,:),z(187+33:189+33,:)))')
ylabel('AC Currents')


%% modify raw measurements

z_new = z;
z_new([190:192,190+33:192+33],:) = z_new([190:192,190+33:192+33],:)/2;
csvwrite('raw_measurements_03-15-2024_15-46.csv',z_new)


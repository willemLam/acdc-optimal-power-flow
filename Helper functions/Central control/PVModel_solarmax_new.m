
function [Pdc] =PVModel_solarmax_new(S, T, nm)
% This function computes the DC power starting from measured GNI and air
% temperature.

% S is the GNI 
% T is air temperature, same dimension of S
% nm is number of modules (51 for the 13 kW Perun and 28 for the 7 kW
% Solarmax)

%% Modules characteristics
Tnoct=45;
Ns=60; %number cell series
Np=1; %number cell parallel
Tref=298.15;
Sref=1000;
Eg_r=1.795842e-19;
q=1.602e-19;    % [C] electronic charge
k=1.3806503e-23;     % [J/K] Boltzmann's constant

Voc_r=37.6;
Isc_r=8.79;

beta= -0.33*Voc_r/100; % [1/K] voltage temp coefficient  change in the solar panel's open-circuit voltage due to an increase in temperature.
alpha= 0.055*Isc_r/100;% [1/K] current T coefficient

%

%% 5 parameters

Iz_r=1.204386402584148e-10;
Rs_r=0.00778;
n_r=0.9514;
Ns=60;
Rp_r=53.07;
Il_r=7.9015;

% 
% Iz_r=1.100e-10;
% Rs_r=0.01; 
% n_r=0.968;
a_ref=k*Tref*Ns*n_r/q;
% Il_r=8.7846;
% Rp_r=48.88;



for i=1:length(S)
    
    
    T(i)=T(i)+273.15;
    
    %Get Tcell from Tair
    T(i)=T(i)+(Tnoct-20).*S(i)/800 + 3*S(i)/800; 

    %Parameters update
    a(i)=a_ref*T(i)/Tref;
    

    Eg(i)=Eg_r*(1-0.0002677*(T(i)-Tref));
    
    Iz(i) = Iz_r*(T(i)/Tref)^3*exp((Eg_r/(k*Tref)-Eg(i)/(k*T(i))));
   
    Il(i)=S(i)/Sref*(Il_r + alpha*(T(i)-Tref));
    
    Rp(i)=Rp_r*Sref/S(i);
    
    Rs(i)=Rs_r;
    
   
    Voc(i)=0.0256*log(S(i))+0.4762 -(-0.000188*log(S(i))+0.003525)*T(i) +Voc_r;
  
  
    
    Vth(i)=k*T(i)/q;
    
    V(:,i)=linspace(0,Voc(i),100);
    
    
    %for each V, from 0 to Voc we can compute I by solving the diode
    %equation
    
    for l = 1:100
        %First way to solve is numerical. It's more precise but with higher
        %computational effort
     %[IDC(l,i)] = fsolve( @(I) -I + Np*Il(i)-Np*Iz(i)*(exp((V(l,i)+I*Ns/Np*Rs(i))/a(i))-1)-(V(l,i)+I*Rs(i)*Ns/Np)/(Ns/Np*Rp(i)), 0);
        
        %Second way is analytical, using Lambert equation (that is however solved numerically).
        %It s proposed by Jain&Kapoor (Solar Energy Materials 2003)
      IDC(l,i)=-V(l,i)/(Rs(i)*Ns) + n_r*Vth(i)/Rs(i)*(-pvl_lambertw(Rs(i)*Iz(i)*Rp(i)*exp((Rp(i)*(Il(i)*Rs(i)+Rs(i)*Iz(i)+V(l,i)/Ns))/(Vth(i)*n_r*(Rs(i)+Rp(i))))) + (Rp(i)*(Rs(i)*Il(i)+Rs(i)*Iz(i)+V(l,i)/Ns))/(n_r*Vth(i)*(Rp(i)+Rs(i))));
        P(l,i)=IDC(l,i)*V(l,i);
    end
    
    [Pmax,Indice] =  max(P); %max for each column
    Pdc(i)=Pmax(i)*nm;
    Vmax=V(Indice,i);
    Imax=IDC(Indice,i);
    
end
%plot(V,IDC)
% ylabel('Current')
% xlabel('Voltage')
% axis([0 inf 0 inf])







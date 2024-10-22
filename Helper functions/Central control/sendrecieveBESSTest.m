% test code to send and recieve from battery using udp
clc;
clear all
addpath(genpath('communication'))

reading_port = 33003;
byte_size = 2048;

[data] = readFromResources(reading_port, byte_size);

% P - kW - negative is charging 
% Q - kVAR - 0
% VSC = 0 for grid following mode, 1 for grid forming mode
%%
P = 20.5;
Q = 0;
VSC = 0;
dest_IP = '192.168.1.73';

dest_port = 34000;

writeToBESS(P, Q, VSC, dest_port, dest_IP)
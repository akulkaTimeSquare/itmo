load('ident_lab2_v05.mat');
a = zad1.a;
b = zad1.b;
w = zad1.w;
Td = 0.1;
set_param('task1_sim/Discrete Transfer Fcn', 'SampleTime', string(Td));
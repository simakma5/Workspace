clear; clc; close all;
%% Initialization
c = 3e8;
mu0 = 4*pi*1e-7;
f0 = 13e9;
H = 0.252e-3;

%% Inductor by a short high-impedance line
% desired value from AWR's iFilter Filter Synthesis wizard
L = 0.45e-9;
% values for a 0.1mm line obtained from TXLine
eps_ef = 1.68602;
ZL = 128.401;

% equation (5.1.44)
LL = c*L/(ZL*sqrt(eps_ef));
% display value (choose whether including via or not)
% disp(['LL = ' num2str(LL) ' m'])
disp(['LL - H = ' num2str(LL-H) ' m'])

%% Capacitors by a short low-impedance lines
% desired value from AWR's iFilter Filter Synthesis wizard
C1 = 0.24e-12;
C2 = 0.44e-12;
% lowest possible impedance for single-mode propagation
fTE10 = f0;
ZC = 2*mu0*H*fTE10;
disp(['ZC = ' num2str(ZC) ' Ohm'])
% value of effective permittivity from TXLines for the given line
eps_ef = 3.425;
% equation from slide 127
LC1 = c*ZC*C1/sqrt(eps_ef);
LC2 = c*ZC*C2/sqrt(eps_ef);
disp(['LC1 = ' num2str(LC1) ' m'])
disp(['LC2 = ' num2str(LC2) ' m'])


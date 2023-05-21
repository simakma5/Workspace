% Array of two out-of-phase dipoles
clear all;
clc;
close all;

d_start = 0.02;
d_end   = 1;
npoints = 100;
d_step  = (d_end-d_start)/npoints;
d       = linspace(d_start,d_end,npoints);

% parameters of two collinear dipoles L1=L2=lambda/2, radius=0
L1    = 0.5;
L2    = L1;
a     = 0;

Z11   = imped(L1,0);
R11   = real(Z11);
X11   = imag(Z11);

for p = 1:npoints;
    Z12(p) = imped(L1,L2,d(p)); % Mutual impedance
end;

R12 = real(Z12);
X12 = imag(Z12);

% driving resistance
Rin = R11 - R12;

% consider input power to the array 100W
% and calculate current at each element
Pin = 100;
I   = sqrt(Pin./(R11-R12));

figure;
hold on;
plot(d,R12,'k');
plot(d,X12,'r');
grid on;
xlabel('d/\lambda'); ylabel('Z_{12} [\Omega]');
legend('R_{12}','X_{12}');
title('Mutual impedance between two side-by-side \lambda/2 dipoles');

figure;
plot(d,Rin,'k');
grid on;
legend('R_{in}');
xlabel('d/\lambda'); ylabel('R_{in} [\Omega]');
title('Input driving resistance');

figure;
plot(d,I,'k');
grid on;
legend('I_{element}');
xlabel('d/\lambda'); ylabel('I_{A} [\Omega]');
title('Current at element for input power to the array 100W');


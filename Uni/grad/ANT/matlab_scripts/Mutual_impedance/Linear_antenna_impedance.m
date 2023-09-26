% Evaluation of input impedance of linear antenna
% 24.2.2021

clear all;
clc;
close all;

L_start = 0.1;
L_end   = 3;
npoints = 100;
L_step  = (L_end-L_start)/npoints;
L       = linspace(L_start,L_end,npoints);

for p = 1:npoints;
        Zin(p)  = imped(L(p),0); % input impedance (at the center)
        Zrad(p) = Zin(p)*(sin(pi*L(p)))^2; % impedance referred to current maximum
end;

figure;
hold on;
plot(L,real(Zin),'k');
plot(L,imag(Zin),'r');
plot(L,real(Zrad),'k--');
plot(L,imag(Zrad),'r--');
grid on;
axis([L_start L_end -700 700]);
xlabel('L/\lambda'); ylabel('Z [\Omega]');
legend('R_{in}','X_{in}','R_r','X_A');
title('Input impedance of infinitely thin dipole');

% Effect of radius on Xin, four different radii in lambda
a1 = 0;
a2 = 1/10000;
a3 = 1/1000;
a4 = 1/100;

L_start = 0.3;
L_end   = 0.7;
npoints = 150;
L_step  = (L_end-L_start)/npoints;
L       = linspace(L_start,L_end,npoints);
for p = 1:npoints;
        Zin_a1(p)  = imped(L(p),a1); % input impedance (at the center)
        Zin_a2(p)  = imped(L(p),a2);
        Zin_a3(p)  = imped(L(p),a3);
        Zin_a4(p)  = imped(L(p),a4);
end;

figure;
hold on;
plot(L,imag(Zin_a1),'k');
plot(L,imag(Zin_a2),'r');
plot(L,imag(Zin_a3),'b');
plot(L,imag(Zin_a4),'m');
line([L_start;L_end],[0;0]) 
grid on;
axis([L_start L_end -500 500]);
xlabel('L/\lambda'); ylabel('X_{in} [\Omega]');
legend('a=0','a=1/10000','a=1/1000','a=1/100');
title('Input reactance of dipoles with different radii');

% generate VSWR and RL plot
Z0     = 75 % impedance of the feeding line
G1     = (Zin_a1-Z0)./(Zin_a1+Z0); % voltage reflection coefficients
G2     = (Zin_a2-Z0)./(Zin_a2+Z0);
G3     = (Zin_a3-Z0)./(Zin_a3+Z0);
G4     = (Zin_a4-Z0)./(Zin_a4+Z0);

S11dB1 = 20.*log10(abs(G1));
S11dB2 = 20.*log10(abs(G2));
S11dB3 = 20.*log10(abs(G3));
S11dB4 = 20.*log10(abs(G4));

VSWR1  = (1+abs(G1))./(1-abs(G1));
VSWR2  = (1+abs(G2))./(1-abs(G2));
VSWR3  = (1+abs(G3))./(1-abs(G3));
VSWR4  = (1+abs(G4))./(1-abs(G4));

figure;
hold on;
plot(L,VSWR1,'k');
plot(L,VSWR2,'r');
plot(L,VSWR3,'b');
plot(L,VSWR4,'m');
axis([L_start L_end 1 5]);
grid on;
xlabel('L/\lambda'); ylabel('VSWR');
legend('a=0','a=1/10000','a=1/1000','a=1/100');
title('VSWR (ref. impedance 75\Omega) of dipoles with different radii');

figure;
hold on;
plot(L,S11dB1,'k');
plot(L,S11dB2,'r');
plot(L,S11dB3,'b');
plot(L,S11dB4,'m');
grid on;
xlabel('L/\lambda'); ylabel('S_{11} in dB');
legend('a=0','a=1/10000','a=1/1000','a=1/100');
title('S_{11} (ref. impedance 75\Omega) of dipoles with different radii');
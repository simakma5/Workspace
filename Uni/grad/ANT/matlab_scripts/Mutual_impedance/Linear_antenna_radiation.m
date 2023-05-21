% Thin-wire antenna of arbitrary length with sinusoidal current
% Evaluate current distribution, radiation pattern plots and directivity
% calculation
% 24.2.2021

clear all;
close all;
clc;

lambda = 1;
k = 2*pi/lambda;
L = (10)*lambda; % length of the dipole

% plot sinusoidal current (I0=1)
z       = linspace(-L/2,L/2,100);
current = sin(k*(L/2-abs(z))); % sin current approximation
%current = 1.*z./z; % constant current
figure;
plot(z./lambda,current);
grid on;
xlabel('z [\lambda]'),ylabel('I [A]');
title('Current distribution along the wire');

% radiation pattern
thdeg  = linspace(1,180,180);
th     = thdeg.*pi/180;

f  = (cos(k*L/2*cos(th))-cos(k*L/2))./sin(th); % field pattern function
fn = f./max(abs(f)); % normalization

% plot pattern in cartesian coordinates
figure
hold on;
plot(thdeg,fn,'b--');
plot(thdeg,abs(fn),'k');
plot(thdeg,(abs(fn)).^2,'r');
line([1;180],[0.5;0.5]);
grid on;
xlabel('\theta [deg]'); ylabel('Pattern level');
legend('f_n','|f_n|','f_n^2');
title('Voltage and power pattern');

% plot in polar coordinates (linear and dB scale)
figure;
abp(th,abs(fn));
figure;
dbp(th,abs(fn),30,20);

% evaluate directivity
integrand = fn.^2.*sin(th);
D    = (4*pi)/(2*pi*trapz(th,integrand))
DdBi = 10*log10(D)


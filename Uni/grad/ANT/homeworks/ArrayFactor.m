%% uniform isotropic array

% broadside array: N=2, d = 1/2, theta_0 = 90
% endfire array: N=2, d = 1/2, theta_0  = 0
% cardidoid: N=2, d = 1/4, theta_0 = 0 / 180

clear all; close all; clc;

theta = linspace(0,360,1360);

theta_0 = 0; % array scan angle in deg (0...endfire, 90..broadside)
d     = 0.4; % array elements separation in lambda
N     = 10;   % number of elements in array
kd    = 2*pi*d; % 2*pi*lambda*d/lambda = 2*pi*d
alpha = -kd*cosd(theta_0);  % progressive phase along the array in rad.
rad2deg(alpha)

% Hansen-Woodyard increased directivity endfire array
alpha = -(kd+pi/N);
rad2deg(alpha)

psi   = kd.*cosd(theta)+alpha;

AF    = sin(N.*psi./2)./(N.*sin(psi./2)); % Array Factor
AFn   = abs(AF./max(AF));

figure;
plot(theta,abs(AFn));
grid on;

figure;
plot(theta,20*log10(abs(AFn)));
grid on;

figure;
theta_radians = deg2rad(theta);
polarplot(theta_radians,abs(AFn));

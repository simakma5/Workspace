% Current distribution and input impedance of the dipole using the method of moments
close all; clear; clc

% number of elements
N = 21;
f = 3e8;
c = 3e8;
lambda = c/f;
w = 2*pi*f;

eps0 = 8.85e-12;

% dipole radius and length
a = 0.005*lambda;   % important: a << lambda
l = 0.47*lambda;    % non-zero thickness influences the ideal lambda/2 model

k=2*pi/lambda;

% element size
h = l/N;

% vector of cell centers
zm = ((1:N)-0.5).*h-(l/2);

% gap voltage
Vo = 1;
mid = (N+1)/2;
V = zeros(N,1);

% gap excitation
V(mid) = -1i*w*eps0*Vo/h;

% Amn matrix elements
A = zeros(N);

for m=1:N
for n=1:N
    % Here, z corresponds to z' (source point).
    r = @(z) sqrt(a^2 + (zm(m)-z).^2);
    f = @(z) exp(-1i*k.*r(z))./(4*pi.*r(z).^5).*((1+1i*k.*r(z)).*(2*r(z).^2-3*a*a)+(k*a.*r(z)).^2);                
    A(m,n) = integral(f,zm(n)-h/2,zm(n)+h/2);
end
end

% Coefficients of current distribution
alpha = A\V;

% input impedance
Zin = Vo/alpha(mid)

figure(1)
plot(abs(alpha))
xlim([0 N+1]);
clear all; close all; clc;
% Two isotropic radiators along X axis - optimize directivity
% array direction set to theta = 90, phi = 0 --> ENDFIRE
p11 = 1;
u11 = p11/(4*pi);
s = linspace(2*pi*0.01,pi,20);

for n = 1:length(s); % s = kd (d=0.01 ... 1 lambda)
    u12 = u11*exp(j*s(n));
    u21 = conj(u12);
    p12 = p11*(sin(s(n))/(s(n)));
    p   = [p11 p12;p12 p11];
    u   = [u11 u12;u21 u11];
    [I, D] = eig(4*pi*u,p)    % eigenvalue equation
    Deig(n)   = D(2,2);
    Ieig = I(:,2);
    Ieig1   = Ieig(1,:);
    Ieig2   = Ieig(2,:);
    PhaseI2I1eig(n) = rad2deg(unwrap(angle(Ieig1/Ieig2))); % relative phase
    
    % solution with inverse of power matrix
    V    = [exp(j*s(n)/2);exp(-j*s(n)/2)];
    Iopt = (1/(4*pi))*inv(p)*V;
    Iopt = [1;-1]; % Force +I, -I currents
    Dopt(n) = 4*pi*(Iopt'*u*Iopt)/(Iopt'*p*Iopt);
    Iopt1   = Iopt(1,:);
    Iopt2   = Iopt(2,:);
    PhaseI2I1(n) = rad2deg(unwrap(angle(Iopt1/Iopt2)));
end;
   
figure;
plot(s/(2*pi),PhaseI2I1eig,'Linewidth',2)
xlabel('s [\lambda]'); ylabel('relative phase of currents');
grid on;

figure;
plot(s/(2*pi),10*log10(Deig),'Linewidth',2)
xlabel('s [\lambda]'); ylabel('D [dBi]');
title('Superdirective currents');
grid on;

figure;
plot(s/(2*pi),10*log10(Dopt),'Linewidth',2)
xlabel('s [\lambda]'); ylabel('D [dBi]');
title('Currents +I, -I');
grid on;

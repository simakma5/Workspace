close all; clear; clc; addpath(genpath(fullfile([pwd '\Uni' '\grad' '\NKA' '\project01'])))

%% 1 Patch antenna TLM (Transmission Line Model) design
% universal constants
c = 3e8;
mu0 = 4*pi*1e-7;
% resonant frequency
f_r = 1.8e9;
lambda0 = c/f_r;
k0 = 2*pi/lambda0;
% copper cladding
t = 35e-6;
sigma = 5.6e7;
% substrate DICLAD
h = 1.5e-3;
epsR = 2.6;
tand = 0.0022;
% computational parameters
fSpan = 1e9;                                    % frequency plots span
Nfpoints = 5e3+1;                               % number of freq. points
f = linspace(f_r-fSpan/2, f_r+fSpan/2, Nfpoints);
NLpoints = 5e3;                                 % number of length points

% 1.1
% antenna dimensions
W = c./(2*f)*sqrt(2/(epsR+1));
epsR_eff = (epsR+1)/2+(epsR-1)/2*(1+12*h./W).^(-1/2);
deltaL = 0.412*h.*(((epsR_eff+0.3).*(W/h+0.264))./((epsR_eff-0.258).*(W/h+0.813)));
L_eff = c./(2*f.*sqrt(epsR_eff));
L = L_eff-2*deltaL;
fprintf('1.1\n----\n')
fprintf('f_r = %.2f GHz\n', f_r*1e-9)
fprintf('W(f=f_r) = %.2f mm\n', W(f==f_r)*1e3)
fprintf('L(f=f_r) = %.2f mm\n', L(f==f_r)*1e3)

% 1.2
% radiating slot admittance Y_s = G + iB
G = W/(120*lambda0)*(1-(k0*h)^2/24);            % for h/lambda0 < 1/10
B = W/(120*lambda0)*(1-0.636*log(k0*h));        % for 0.35 < W/lambda0 < 2
Y_s = G + 1i*B;
Z_s = 1./Y_s;

% characteristic impedance Z_c of a microstrip
if W/h <= 1
    Z_c = 60/sqrt(epsR_eff)*log(8*h/W+W/(4*h));
else
    Z_c = 120*pi./(sqrt(epsR_eff).*(W/h+1.393+0.667*log(W/h+1.444)));
end

% input impedance Z_in as a function of:
% L ... position of the feeding point from the edge of the antenna,
% f ... operating frequency
beta = 2*pi*sqrt(epsR_eff)/lambda0;
[L1, L2, Z1, Z2] = deal(zeros(NLpoints, Nfpoints));
for fi = 1:Nfpoints
    L1(:,fi) = linspace(0, L(fi), NLpoints);
    L2(:,fi) = L(fi) - L1(:, fi);
    Z1(:,fi) = Z_c(fi).*(Z_s(fi)+1i*Z_c(fi).*tan(beta(fi).*L1(:,fi)))./(Z_c(fi)+1i.*Z_s(fi).*tan(beta(fi).*L1(:,fi)));
    Z2(:,fi) = Z_c(fi).*(Z_s(fi)+1i*Z_c(fi).*tan(beta(fi).*L2(:,fi)))./(Z_c(fi)+1i.*Z_s(fi).*tan(beta(fi).*L2(:,fi)));
end
Z_in = 1./(1./Z1+1./Z2);

% input impedance frequency plot for L1 == 0 mm
figure(1)
plot(f*1e-9, real(Z_in(1,:)))
hold on
plot(f*1e-9, imag(Z_in(1,:)))
xline(f_r*1e-9, '-.')
yline(0)
hold off
axis tight
grid on
grid minor
xlabel('Frequency [GHz]')
ylabel('[\Omega]')
legend('R_{in}', 'X_{in}')
title('Input impedance of edge-fed antenna')

% module of reflection coefficient frequency plot for L1 == 0 mm
figure(2)
plot(f*1e-9, 20*log10(abs((50-Z_in(1,:))./(50+Z_in(1,:)))))
hold on
xline(f_r*1e-9, '-.')
hold off
axis tight
grid on
grid minor
xlabel('Frequency [GHz]')
ylabel('|\Gamma| [dB]')
title('Reflection coefficient of edge-fed antenna')

% 1.3
% input impedance for L1 == 0 mm, f == f_r
fprintf('\n1.3\n----\n')
fprintf('Z_in(L1=0,f=f_r) = %.2f%+.2fi Ohm\n', real(Z_in(1,f==f_r)), imag(Z_in(1,f==f_r)))

% input impedance feeding point position plot
figure(3)
plot(L1(:,f==f_r)/L(f==f_r), real(Z_in(:,f==f_r)))
hold on
plot(L1(:,f==f_r)/L(f==f_r), imag(Z_in(:,f==f_r)))
line = yline(50, '--', '50');
line.LabelVerticalAlignment = 'bottom';
hold off
axis tight
grid on
grid minor
xlabel('L_1/L [-]')
ylabel('[\Omega]')
legend('R_{in}', 'X_{in}', 'Location', 'southeast')
title('Input impedance at resonant frequency')

% find L for which R_in is approx. 50 Ohm @f==f_r
[R_match, l_match] = min(abs(real(Z_in(:,f==f_r))-50));
R_match = R_match + 50;
disp(['L_match = ' num2str(round(L1(l_match,f==f_r)*1e3, 2)) ' mm'])
% disp(['R_match = ' num2str(round(R_match, 2)) ' Ohm'])

% input impedance frequency plot for L == L_match
figure(4)
plot(f*1e-9, real(Z_in(l_match,:)))
hold on
plot(f*1e-9, imag(Z_in(l_match,:)))
xline(f_r*1e-9, '-.')
hold off
axis tight
grid on
grid minor
xlabel('Frequency [GHz]')
ylabel('[\Omega]')
legend('R_{in}', 'X_{in}')
title('Input impedance of matched antenna')

% module of reflection coefficient frequency plot for L == L_match
figure(5)
Gamma_match_dB = 20*log10(abs((50-Z_in(l_match,:))./(50+Z_in(l_match,:))));
plot(f*1e-9, Gamma_match_dB)
hold on
xline(f_r*1e-9, '-.')
hold off
axis tight
grid on
grid minor
xlabel('Frequency [GHz]')
ylabel('|\Gamma| [dB]')
title('Reflection coefficient of matched antenna')

% 1.4
% Q-factor defined as the ratio of energy accumulated in cavity over
% energy lost in materials
Q_d = 1/tand;                                   % dielectric
Q_c = sqrt(pi*f_r*mu0*sigma*t);                 % conductor
Q_r = c*sqrt(epsR)/(4*f_r*h);                   % radiation
Q_t = 1/(1/Q_d+1/Q_c+1/Q_r);                    % total

% analytical formula for bandwidth of VSWR < 2
% BW = (VSWR-1)/(Q_t*sqrt(VSWR))
BW_analytic = 1/(Q_t*sqrt(2))*100;
fprintf('\n1.4\n----\n')
fprintf('Analytically: BW(VSWR<2) = %.2f %%\n', BW_analytic)

% determine BW(VSWR<2) (|Gamma(VSWR=2)| = -9.54 dB) from computed |Gamma|
[~, f_VSWR2] = min(abs(Gamma_match_dB+9.54));
BW_numeric = abs(Nfpoints-2*f_VSWR2)/Nfpoints*100;
fprintf('Numerically:  BW(VSWR<2) = %.2f %%\n', BW_numeric)

% 1.5
% evaluate BW(VSWR<2) according to the analytical formula for commonly
% manufactured ranges of values of h and epsR
hRange = 0.5e-3:0.1e-3:5e-3;
epsRRange = 1:0.1:10;
Q_d = 1/tand;
Q_c = sqrt(pi*f_r*mu0*sigma*t);
Q_r = zeros(length(epsRRange), length(hRange));
for i = 1:length(epsRRange)
    Q_r(i,:) = c*sqrt(epsRRange(i))./(4*f_r*hRange(:));
end
Q_t = 1./(1/Q_d+1/Q_c+1./Q_r);
BW_analytic = 1./(Q_t*sqrt(2))*100;

figure(6)
plot(hRange*1e3, BW_analytic(epsRRange==epsR,:))
line = xline(h*1e3, '--', 'h');
line.LabelVerticalAlignment = 'bottom';
axis tight
grid on
grid minor
xlabel('h [mm]')
ylabel('BW(VSWR<2) [%]')
title(['Antenna bandwidth of VSWR < 2 for epsR = ' num2str(epsR)])

figure(7)
plot(epsRRange, BW_analytic(:,hRange==h))
line = xline(epsR, '--', 'epsR');
line.LabelVerticalAlignment = 'bottom';
axis tight
grid on
grid minor
xlabel('epsR [-]')
ylabel('BW(VSWR<2) [%]')
title(['Antenna bandwidth of VSWR < 2 for h = ' num2str(h*1e3) ' mm'])

% 1.6
theta = linspace(-pi/2, pi/2, 180);
F_E = sin(k0*h/2.*cos(theta))./(k0*h/2.*cos(theta)).*cos(k0*L_eff(f==f_r)/2.*sin(theta));
F_H = cos(theta).*sin(k0*h/2.*cos(theta))./(k0*h/2.*cos(theta)).*sin(k0*W(f==f_r)/2.*sin(theta))./(k0*W(f==f_r)/2.*sin(theta));

figure(81)
polaraxes( ...
    'RLim', [-40 0], ...
    'ThetaLim', [-90 90], ...
    'thetadir', 'clockwise', ...
    'thetazerolocation', 'top');
hold on
polarplot(theta, 10*log10(F_E)/max(F_E))
hold off
title('Normalized elevation field pattern')

figure(91)
polaraxes( ...
    'RLim', [-40 0], ...
    'ThetaLim', [-90 90], ...
    'thetadir', 'clockwise', ...
    'thetazerolocation', 'top');
hold on
polarplot(theta, 10*log10(F_H)/max(F_H))
hold off
title('Normalized azimuthal field pattern')

%% 2 Patch antenna measurement
patch4mm = sparameters(fullfile('data', '01 PATCH.S1P'));
patch16mm = sparameters(fullfile('data', '03 PATCH.S1P'));
patch28mm = sparameters(fullfile('data', '05 PATCH.S1P'));
% fMeas = patch4mm.Frequencies;
% NfpointsMeas = length(fMeas);

L4mm = find(round(L1(:,f==f_r),3)==4e-3);
L4mm = L4mm(ceil(end/2));
L16mm = find(round(L1(:,f==f_r),3)==16e-3);
L16mm = L16mm(ceil(end/2));
L26mm = find(round(L1(:,f==f_r),3)==26e-3);
L26mm = L26mm(ceil(end/2));

figure(11)
rfplot(patch4mm)
% rfplot(patch16mm)
% rfplot(patch28mm)
hold on
plot(f*1e-9, 20*log10(abs((50-Z_in(L4mm,:))./(50+Z_in(L4mm,:)))))
hold off
xlabel('Frequency [GHz]')
ylabel('|S_{11}| [dB]')
legend('Measured data', 'TLM', 'Location', 'southeast')
title('')

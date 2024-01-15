close all; clear; clc; addpath(genpath(fullfile([pwd '\Uni' '\grad' '\NKA' '\project1'])))
% s = settings;
% s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";
% universal constants
c = 299792458;
mu0 = 4*pi*1e-7;
% resonant frequency
f_r = 1.8e9;
lambda0 = c/f_r;
k0 = 2*pi/lambda0;

%% 1 Patch antenna TLM (Transmission Line Model) design
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
NLpoints = 1e3;                                 % number of length points

% 1.1
% antenna dimensions
W = c/(2*f_r)*sqrt(2/(epsR+1));
epsR_eff = (epsR+1)/2+(epsR-1)/2*(1+12*h/W)^(-1/2);
deltaL = 0.412*h*(((epsR_eff+0.3)*(W/h+0.264))/((epsR_eff-0.258)*(W/h+0.813)));
L_eff = c/(2*f_r*sqrt(epsR_eff));
L = L_eff-2*deltaL;
fprintf('1.1\n----\n')
fprintf('f_r = %.2f GHz\n', f_r*1e-9)
fprintf('W(f=f_r) = %.2f mm\n', W*1e3)
fprintf('L(f=f_r) = %.2f mm\n', L*1e3)

% 1.2
% input impedance frequency plot for feed offset of 0 mm
Z_in = impedanceTLM(epsR, h, W, L, 0, f);
figure('Name', 'Input impedance of an edge-fed antenna')
plot(f*1e-9, real(Z_in))
hold on
plot(f*1e-9, imag(Z_in))
xline(f_r*1e-9, '-.')
yline(0)
hold off
axis tight
grid on
grid minor
xlabel('Frequency [GHz]')
ylabel('[\Omega]')
legend('R_{in}', 'X_{in}')
% title('Input impedance of an edge-fed antenna')

% module of reflection coefficient frequency plot for L1 == 0 mm
figure('Name', 'Reflection coefficient of an edge-fed antenna')
plot(f*1e-9, 20*log10(abs((50-Z_in)./(50+Z_in))))
hold on
xline(f_r*1e-9, '-.')
hold off
axis tight
grid on
grid minor
xlabel('Frequency [GHz]')
ylabel('|\Gamma| [dB]')
% title('Reflection coefficient of an edge-fed antenna')

% 1.3
% input impedance for L1 == 0 mm && f == f_r
fprintf('\n1.3\n----\n')
fprintf('Z_in(L1=0,f=f_r) = %.2f%+.2fi Ohm\n', real(Z_in(f==f_r)), imag(Z_in(f==f_r)))

% input impedance feeding point position plot for f == f_r
L1 = linspace(0, L, NLpoints);
Z_in = impedanceTLM(epsR, h, W, L, L1, f_r);
figure('Name', 'Input impedance at resonance of an offset-fed antenna')
plot(L1/L, real(Z_in))
hold on
plot(L1/L, imag(Z_in))
line = yline(50, '--', '50');
line.LabelVerticalAlignment = 'bottom';
hold off
axis tight
grid on
grid minor
xlabel('L_1/L [-]')
ylabel('[\Omega]')
legend('R_{in}', 'X_{in}', 'Location', 'southeast')
% title('Input impedance at resonant frequency')

% find L for which R_in is approx. 50 Ohm at f == f_r
[R_match, L1_match] = min(abs(real(Z_in)-50));
R_match = R_match + 50;
disp(['L_match = ' num2str(round(L1(L1_match)*1e3, 2)) ' mm'])
% disp(['R_match = ' num2str(round(R_match, 2)) ' Ohm'])

% input impedance frequency plot for L == L_match
Z_in = impedanceTLM(epsR, h, W, L, L1(L1_match), f);
figure('Name', 'Input impedance of a matched antenna')
plot(f*1e-9, real(Z_in))
hold on
plot(f*1e-9, imag(Z_in))
xline(f_r*1e-9, '-.')
hold off
axis tight
grid on
grid minor
xlabel('Frequency [GHz]')
ylabel('[\Omega]')
legend('R_{in}', 'X_{in}')
% title('Input impedance of a matched antenna')

% module of reflection coefficient frequency plot for L == L_match
figure('Name', 'Reflection coefficient of a matched antenna')
Gamma_match_dB = 20*log10(abs((50-Z_in./(50+Z_in))));
plot(f*1e-9, Gamma_match_dB)
hold on
xline(f_r*1e-9, '-.')
hold off
axis tight
grid on
grid minor
xlabel('Frequency [GHz]')
ylabel('|\Gamma| [dB]')
% title('Reflection coefficient of a matched antenna')

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

figure('Name', ['Antenna bandwidth (VSWR < 2) for epsR = ' num2str(epsR)])
plot(hRange*1e3, BW_analytic(epsRRange==epsR,:))
line = xline(h*1e3, '--', 'h');
line.LabelVerticalAlignment = 'bottom';
axis tight
grid on
grid minor
xlabel('h [mm]')
ylabel('BW(VSWR<2) [%]')
% title(['Antenna bandwidth of VSWR < 2 for epsR = ' num2str(epsR)])

figure('Name', ['Antenna bandwidth (VSWR < 2) for h = ' num2str(h*1e3) ' mm'])
plot(epsRRange, BW_analytic(:,hRange==h))
line = xline(epsR, '--', 'epsR');
line.LabelVerticalAlignment = 'bottom';
axis tight
grid on
grid minor
xlabel('epsR [-]')
ylabel('BW(VSWR<2) [%]')
% title(['Antenna bandwidth of VSWR < 2 for h = ' num2str(h*1e3) ' mm'])

% 1.6
theta = linspace(-pi/2, pi/2, 180);
F_E = sin(k0*h/2.*cos(theta))./(k0*h/2.*cos(theta)).*cos(k0*L_eff/2.*sin(theta));
F_H = cos(theta).*sin(k0*h/2.*cos(theta))./(k0*h/2.*cos(theta)).*sin(k0*W/2.*sin(theta))./(k0*W/2.*sin(theta));

figure('Name', 'Normalized elevation field pattern')
polaraxes( ...
    'RLim', [-40 0], ...
    'ThetaLim', [-90 90], ...
    'thetadir', 'clockwise', ...
    'thetazerolocation', 'top');
hold on
polarplot(theta, 10*log10(F_E)/max(F_E))
hold off
% title('Normalized elevation field pattern')

figure('Name', 'Normalized azimuthal field pattern')
polaraxes( ...
    'RLim', [-40 0], ...
    'ThetaLim', [-90 90], ...
    'thetadir', 'clockwise', ...
    'thetazerolocation', 'top');
hold on
polarplot(theta, 10*log10(F_H)/max(F_H))
hold off
% title('Normalized azimuthal field pattern')

%% 2 Patch antenna measurement
% % load measurement S-parameter data
% S11_4mm = sparameters(fullfile('data', '01 PATCH.S1P'));
% S11_16mm = sparameters(fullfile('data', '03 PATCH.S1P'));
% S11_28mm = sparameters(fullfile('data', '05 PATCH.S1P'));

% measurement parameters
f_meas = sparameters(fullfile('data', '01 PATCH.S1P')).Frequencies;
lambda0_meas = c./f_meas;

% measured antenna properties
epsR_meas = 1;
h_meas = 5e-3;
L_meas = 58e-3;
W_meas = 61e-3;

L1_meas = [4e-3; 16e-3; 28e-3];
S11_meas = [
    squeeze(sparameters(fullfile('data', '01 PATCH.S1P')).Parameters(1,1,:)), ...
    squeeze(sparameters(fullfile('data', '03 PATCH.S1P')).Parameters(1,1,:)), ...
    squeeze(sparameters(fullfile('data', '05 PATCH.S1P')).Parameters(1,1,:))
];
Z_in_meas = [
    impedanceTLM(epsR_meas, h_meas, W_meas, L_meas, L1_meas(1), S11_meas(1)), ...
    impedanceTLM(epsR_meas, h_meas, W_meas, L_meas, L1_meas(2), S11_meas(2)), ...
    impedanceTLM(epsR_meas, h_meas, W_meas, L_meas, L1_meas(3), S11_meas(3))
];

% continue from editing the following line and on
z_meas = ((1 + s11(:, 2:4))./(1 - s11(:, 2:4)))*z0;

% figure(11)
% rfplot(S11_4mm)
% hold on
% plot(f*1e-9, 20*log10(abs((50-Z_in(L4mm,:))./(50+Z_in(L4mm,:)))))
% hold off
% xlabel('Frequency [GHz]')
% ylabel('|S_{11}| [dB]')
% legend('Measured data', 'TLM', 'Location', 'southeast')
% title('')

%% Functions
% Calculation of input impedance using the Transmission Line Model
function Z_in = impedanceTLM(epsR, h, W, L, feedOffset, f)
    % universal quantities
    c = 299792458;
    lambda0 = c./f;
    k0 = 2*pi./lambda0;

    % effective permittivity
    epsR_eff = (epsR+1)/2+(epsR-1)/2*(1+12*h./W).^(-1/2);

    % radiating slot admittance Y_s = G + iB
    G = W./(120*lambda0).*(1-(k0*h).^2/24);             % for h/lambda0 < 1/10
    B = W./(120*lambda0).*(1-0.636*log(k0*h));          % for 0.35 < W/lambda0 < 2
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
    beta = 2*pi*sqrt(epsR_eff)./lambda0;
    Z1 = Z_c.*(Z_s+1i*Z_c.*tan(beta*feedOffset))./(Z_c+1i*Z_s.*tan(beta*feedOffset));
    Z2 = Z_c.*(Z_s+1i*Z_c.*tan(beta*(L-feedOffset)))./(Z_c+1i*Z_s.*tan(beta*(L-feedOffset)));
    Z_in = 1./(1./Z1+1./Z2);
end
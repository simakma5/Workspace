close all; clear; clc; addpath(genpath(fullfile([pwd, '\Uni', '\grad', '\NKA', '\project1', '\script'])))
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";
% universal constants
c = 299792458;
mu0 = 4*pi*1e-7;
% reference impedance
Z0 = 50;

%% 1 Patch antenna TLM (Transmission Line Model) design
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

%% 1.1 Antenna dimensions
W = c/(2*f_r)*sqrt(2/(epsR+1));
epsR_eff = (epsR+1)/2+(epsR-1)/2*(1+12*h/W)^(-1/2);
deltaL = 0.412*h*(((epsR_eff+0.3)*(W/h+0.264))/((epsR_eff-0.258)*(W/h+0.813)));
L_eff = c/(2*f_r*sqrt(epsR_eff));
L = L_eff-2*deltaL;
fprintf('1.1 Antenna dimensions\n----\n')
fprintf('f_r = %.2f GHz\n', f_r*1e-9)
fprintf('W(f=f_r) = %.2f mm\n', W*1e3)
fprintf('L(f=f_r) = %.2f mm\n', L*1e3)

%% 1.2 Input impedance and reflection coefficient frequency dependence
% input impedance and reflection coefficient frequency plot for feed offset of 0 mm
Z_in = impedanceTLM(epsR, h, W, L, 0, f);
figure('Name', 'Edge-fed antenna (zero feed offset) input impedance and reflection coefficient', 'Units', 'normalized', 'Position', [0, 0, 1, 1])
tiledlayout(1,2)
nexttile;
    plot(f*1e-9, real(Z_in))
    hold on
    plot(f*1e-9, imag(Z_in))
    xline(f_r*1e-9, '-.')
    yline(0)
    hold off
    grid on
    grid minor
    xlim([f(1)*1e-9, f(end)*1e-9])
    xlabel('Frequency [GHz]')
    ylabel('Z_{in} [\Omega]')
    legend('R_{in}', 'X_{in}')
nexttile;
    plot(f*1e-9, 20*log10(abs((Z0-Z_in)./(Z0+Z_in))))
    hold on
    xline(f_r*1e-9, '-.')
    hold off
    grid on
    grid minor
    xlim([f(1)*1e-9, f(end)*1e-9])
    ylim([-30, 0])
    xlabel('Frequency [GHz]')
    ylabel('|\Gamma| [dB]')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 26)
saveas(gcf, fullfile([pwd, '\Uni', '\grad', '\NKA', '\project1', '\report', '\src', '\edge-fed-antenna']), 'epsc')

%% 1.3 Feed-offset antenna matching
% input impedance for L1 == 0 mm && f == f_r
fprintf('\n1.3 Feed-offset antenna matching\n----\n')
fprintf('Z_in(L1=0,f=f_r) = %.2f%+.2fi Ohm\n', real(Z_in(f==f_r)), imag(Z_in(f==f_r)))

% input impedance and reflection coefficient feed-offset plot for f == f_r
L1 = linspace(0, L, NLpoints);
Z_in = impedanceTLM(epsR, h, W, L, L1, f_r);
figure('Name', 'Input impedance at resonance of an offset-fed antenna', 'Units', 'normalized', 'Position', [0, 0, 1, 1])
plot(L1/L, real(Z_in))
hold on
plot(L1/L, imag(Z_in))
line = yline(Z0, '--', '50');
line.LabelVerticalAlignment = 'bottom';
hold off
grid on
grid minor
xlabel('L_1/L [-]')
ylabel('Z_{in} [\Omega]')
legend('R_{in}', 'X_{in}', 'Location', 'southeast')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 26)
saveas(gcf, fullfile([pwd, '\Uni', '\grad', '\NKA', '\project1', '\report', '\src', '\feed-offset-impedance-plot']), 'epsc')

% find L for which R_in is approx. 50 Ohm at f == f_r
[R_match, L1_matchIdx] = min(abs(real(Z_in)-Z0));
R_match = R_match + Z0;
disp(['L_match = ' num2str(round(L1(L1_matchIdx)*1e3, 2)) ' mm'])
disp(['R_match = ' num2str(round(R_match, 2)) ' Ohm'])

% input impedance and reflection coefficient frequency plot for feed offset of 0 mm
Z_in = impedanceTLM(epsR, h, W, L, L1(L1_matchIdx), f);
fprintf('Z_in(L1=L_match,f=f_r) = %.2f%+.2fi Ohm\n', real(Z_in(f==f_r)), imag(Z_in(f==f_r)))
figure('Name', 'Matched antenna input impedance and reflection coefficient', 'Units', 'normalized', 'Position', [0, 0, 1, 1])
tiledlayout(1,2)
nexttile;
    plot(f*1e-9, real(Z_in))
    hold on
    plot(f*1e-9, imag(Z_in))
    xline(f_r*1e-9, '-.')
    hold off
    grid on
    grid minor
    xlim([f(1)*1e-9, f(end)*1e-9])
    xlabel('Frequency [GHz]')
    ylabel('Z_{in} [\Omega]')
    legend('R_{in}', 'X_{in}')
nexttile;
    plot(f*1e-9, 20*log10(abs((Z0-Z_in)./(Z0+Z_in))))
    hold on
    xline(f_r*1e-9, '-.')
    hold off
    grid on
    grid minor
    xlim([f(1)*1e-9, f(end)*1e-9])
    ylim([-30, 0])
    xlabel('Frequency [GHz]')
    ylabel('|\Gamma| [dB]')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 26)
saveas(gcf, fullfile([pwd, '\Uni', '\grad', '\NKA', '\project1', '\report', '\src', '\matched-antenna']), 'epsc')

%% 1.4 Antenna bandwidth of VSWR < 2
% Q-factor defined as the ratio of energy accumulated in cavity over
% energy lost in materials
Q_d = 1/tand;                                   % dielectric
Q_c = sqrt(pi*f_r*mu0*sigma*t);                 % conductor
Q_r = c*sqrt(epsR)/(4*f_r*h);                   % radiation
Q_t = 1/(1/Q_d+1/Q_c+1/Q_r);                    % total

% analytical formula for bandwidth of VSWR < 2
% BW = (VSWR-1)/(Q_t*sqrt(VSWR))
BW_analytic = 1/(Q_t*sqrt(2))*100;
fprintf('\n1.4 Antenna bandwidth of VSWR < 2\n----\n')
fprintf('Analytically: BW(VSWR<2) = %.2f %%\n', BW_analytic)

% determine BW(VSWR<2) (|Gamma(VSWR=2)| = -9.54 dB) from computed Gamma
[~, fIdx] = min(abs(20*log10(abs((Z0-Z_in./(Z0+Z_in))))+9.54));
BW_numeric = abs(f_r-f(fIdx))/(f_r+f(fIdx))*100;
fprintf('Numerically:  BW(VSWR<2) = %.2f %%\n', BW_numeric)

%% 1.5 Analytical BW of VSWR < 2 for different substrates
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

figure('Name', 'Antenna bandwidth of VSWR < 2', 'Units', 'normalized', 'Position', [0, 0, 1, 1])
tiledlayout(1,2)
nexttile;
    plot(hRange*1e3, BW_analytic(epsRRange==epsR,:))
    line = xline(h*1e3, '--', 'h');
    line.LabelVerticalAlignment = 'bottom';
    grid on
    grid minor
    axis tight
    xlabel('h [mm]')
    ylabel('BW(VSWR<2) [%]')
    title(['epsR = ' num2str(epsR)])
nexttile;
    plot(epsRRange, BW_analytic(:,hRange==h))
    line = xline(epsR, '--', 'epsR');
    line.LabelVerticalAlignment = 'bottom';
    grid on
    grid minor
    axis tight
    xlabel('epsR [-]')
    ylabel('BW(VSWR<2) [%]')
    title(['h = ' num2str(h*1e3) ' mm'])
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 26)
saveas(gcf, fullfile([pwd, '\Uni', '\grad', '\NKA', '\project1', '\report', '\src', '\bandwidth-material-dependence']), 'epsc')

%% 1.6 Antenna field patterns
theta = linspace(-pi/2, pi/2, 180);
F_E = sin(k0*h/2.*cos(theta))./(k0*h/2.*cos(theta)).*cos(k0*L_eff/2.*sin(theta));
F_H = cos(theta).*sin(k0*h/2.*cos(theta))./(k0*h/2.*cos(theta)).*sin(k0*W/2.*sin(theta))./(k0*W/2.*sin(theta));

figure('Name', 'Normalized field patterns', 'Units', 'normalized', 'Position', [0, 0, 1, 1])
tiles = tiledlayout(1,2);
    pax = polaraxes(tiles);
    pax.RLim = [-40, 0];
    pax.ThetaLim = [-90, 90];
    pax.ThetaDir = 'clockwise';
    pax.ThetaZeroLocation = 'top';
    pax.RTickLabels = {'-40 dB', '-30 dB','-20 dB','-10 dB','0 dB'};
    pax.Layout.Tile = 1;
    hold on
    polarplot(theta, 10*log10(F_E)/max(F_E))
    hold off
    title('E-plane (azimuthal)')

    pax = polaraxes(tiles);
    pax.RLim = [-40, 0];
    pax.ThetaLim = [-90, 90];
    pax.ThetaDir = 'clockwise';
    pax.ThetaZeroLocation = 'top';
    pax.RTickLabels = {'-40 dB', '-30 dB','-20 dB','-10 dB','0 dB'};
    pax.Layout.Tile = 2;
    hold on
    polarplot(theta, 10*log10(F_H)/max(F_H))
    hold off
    title('H-plane (elevation)')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 26)
saveas(gcf, fullfile([pwd, '\Uni', '\grad', '\NKA', '\project1', '\report', '\src', '\field-patterns']), 'epsc')

%% 2 Patch antenna measurement
% measurement parameters
f_meas = sparameters(fullfile('data', '01 PATCH.S1P')).Frequencies;
lambda0_meas = c./f_meas;
epsR_meas = 1;
tand = 0;

% measured antenna dimensions
h_meas = 5e-3;
L_meas = 58e-3;
W_meas = 61e-3;

%% Measured vs TLM input impedance and reflection coefficient
% measurement #1
L1_meas = 4e-3;
Gamma_meas = squeeze(sparameters(fullfile('data', '01 PATCH.S1P')).Parameters(1,1,:));
Z_meas = Z0*(1+Gamma_meas)./(1-Gamma_meas);

% model
Z_TLM = impedanceTLM(epsR_meas, h_meas, W_meas, L_meas, L1_meas, f_meas);
Gamma_TLM = (Z_TLM-Z0)./(Z_TLM+Z0);

% plot
figure('Name', 'Measurement 1 input impedance and reflection coefficient', 'Units', 'normalized', 'Position', [0, 0, 1, 1])
tiledlayout(1,2)
nexttile;
    plot(f_meas*1e-9, real(Z_meas))
    hold on
    plot(f_meas*1e-9, imag(Z_meas))
    plot(f_meas*1e-9, real(Z_TLM))
    plot(f_meas*1e-9, imag(Z_TLM))
    yline([50 0], '--')
    hold off
    grid on
    grid minor
    xlabel('Frequency [GHz]')
    ylabel('Z_{in} [\Omega]')
    legend('R_{in} (meas.)', 'X_{in} (meas.)', 'R_{in} (TLM)', 'X_{in} (TLM)', 'location', 'best')
nexttile;
    plot(f_meas*1e-9, 20*log10(abs(Gamma_meas)))
    hold on
    plot(f_meas*1e-9, 20*log10(abs(Gamma_TLM)))
    hold off
    grid on
    grid minor
    xlabel('Frequency [GHz]')
    ylabel('|\Gamma| [dB]')
    legend('Measurement', 'TLM', 'location', 'best')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 26)
saveas(gcf, fullfile([pwd, '\Uni', '\grad', '\NKA', '\project1', '\report', '\src', '\measurement1']), 'epsc')

% measurement #2 (best match)
L1_meas = 16e-3;
Gamma_meas = squeeze(sparameters(fullfile('data', '03 PATCH.S1P')).Parameters(1,1,:));
Z_meas = Z0*(1+Gamma_meas)./(1-Gamma_meas);

% model
Z_TLM = impedanceTLM(epsR_meas, h_meas, W_meas, L_meas, L1_meas, f_meas);
Gamma_TLM = (Z_TLM-Z0)./(Z_TLM+Z0);

% plot
figure('Name', 'Measurement 2 input impedance and reflection coefficient', 'Units', 'normalized', 'Position', [0, 0, 1, 1])
tiledlayout(1,2)
nexttile;
    plot(f_meas*1e-9, real(Z_meas))
    hold on
    plot(f_meas*1e-9, imag(Z_meas))
    plot(f_meas*1e-9, real(Z_TLM))
    plot(f_meas*1e-9, imag(Z_TLM))
    yline([50 0], '--')
    hold off
    grid on
    grid minor
    xlabel('Frequency [GHz]')
    ylabel('Z_{in} [\Omega]')
    legend('R_{in} (meas.)', 'X_{in} (meas.)', 'R_{in} (TLM)', 'X_{in} (TLM)', 'location', 'best')
nexttile;
    plot(f_meas*1e-9, 20*log10(abs(Gamma_meas)))
    hold on
    plot(f_meas*1e-9, 20*log10(abs(Gamma_TLM)))
    hold off
    grid on
    grid minor
    xlabel('Frequency [GHz]')
    ylabel('|\Gamma| [dB]')
    legend('Measurement', 'TLM', 'location', 'best')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 26)
saveas(gcf, fullfile([pwd, '\Uni', '\grad', '\NKA', '\project1', '\report', '\src', '\measurement2']), 'epsc')

% measurement #3
L1_meas = 28e-3;
Gamma_meas = squeeze(sparameters(fullfile('data', '05 PATCH.S1P')).Parameters(1,1,:));
Z_meas = Z0*(1+Gamma_meas)./(1-Gamma_meas);

% model
Z_TLM = impedanceTLM(epsR_meas, h_meas, W_meas, L_meas, L1_meas, f_meas);
Gamma_TLM = (Z_TLM-Z0)./(Z_TLM+Z0);

% plot
figure('Name', 'Measurement 3 input impedance and reflection coefficient', 'Units', 'normalized', 'Position', [0, 0, 1, 1])
tiledlayout(1,2)
nexttile;
    plot(f_meas*1e-9, real(Z_meas))
    hold on
    plot(f_meas*1e-9, imag(Z_meas))
    plot(f_meas*1e-9, real(Z_TLM))
    plot(f_meas*1e-9, imag(Z_TLM))
    yline([50 0], '--')
    hold off
    grid on
    grid minor
    xlabel('Frequency [GHz]')
    ylabel('Z_{in} [\Omega]')
    legend('R_{in} (meas.)', 'X_{in} (meas.)', 'R_{in} (TLM)', 'X_{in} (TLM)', 'location', 'best')
nexttile;
    plot(f_meas*1e-9, 20*log10(abs(Gamma_meas)))
    hold on
    plot(f_meas*1e-9, 20*log10(abs(Gamma_TLM)))
    hold off
    grid on
    grid minor
    xlabel('Frequency [GHz]')
    ylabel('|\Gamma| [dB]')
    legend('Measurement', 'TLM', 'location', 'best')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 26)
saveas(gcf, fullfile([pwd, '\Uni', '\grad', '\NKA', '\project1', '\report', '\src', '\measurement3']), 'epsc')

%% 2.2 Measured and TLM resonant frequencies and BW of VSWR < 2
% find resonant frequency from the measurement data
Gamma = squeeze(sparameters(fullfile('data', '03 PATCH.S1P')).Parameters(1,1,:));
[~, fIdx] = min(abs(Gamma));
f_r = f_meas(fIdx);
fprintf('\n2.2 Measured and TLM resonant frequencies and BW of VSWR < 2\n----\n')
fprintf('(Measurement)\n')
fprintf('Resonant frequency: %.2f\n', f_r*1e-9)

Q_d = 1/tand;                                   % dielectric
Q_c = sqrt(pi*f_r*mu0*sigma*t);            % conductor
Q_r = c*sqrt(epsR)/(4*f_r*h);              % radiation
Q_t = 1/(1/Q_d+1/Q_c+1/Q_r);                    % total

% analytical formula for bandwidth of VSWR < 2
% BW = (VSWR-1)/(Q_t*sqrt(VSWR))
BW_analytic = 1/(Q_t*sqrt(2))*100;
fprintf('Analytically: BW(VSWR<2) = %.2f %%\n', BW_analytic)

% determine BW(VSWR<2) (|Gamma(VSWR=2)| = -9.54 dB) from computed Gamma
[~, fIdx] = min(abs(20*log10(abs(Gamma))+9.54));
BW_numeric = abs(f_r-f_meas(fIdx))/(f_r+f_meas(fIdx))*100;
fprintf('Numerically:  BW(VSWR<2) = %.2f %%\n', BW_numeric)

% find TLM resonant frequency
L1_meas = 16e-3;
Z = impedanceTLM(epsR_meas, h_meas, W_meas, L_meas, L1_meas, f_meas);
Gamma = (Z-Z0)./(Z+Z0);
[~, fIdx] = min(abs(Gamma));
f_r = f_meas(fIdx);
fprintf('\n(TLM)\n')
fprintf('Resonant frequency: %.2f\n', f_r*1e-9)

Q_d = 1/tand;                                   % dielectric
Q_c = sqrt(pi*f_r*mu0*sigma*t);            % conductor
Q_r = c*sqrt(epsR)/(4*f_r*h);              % radiation
Q_t = 1/(1/Q_d+1/Q_c+1/Q_r);                    % total

% analytical formula for bandwidth of VSWR < 2
% BW = (VSWR-1)/(Q_t*sqrt(VSWR))
BW_analytic = 1/(Q_t*sqrt(2))*100;
fprintf('Analytically: BW(VSWR<2) = %.2f %%\n', BW_analytic)

% determine BW(VSWR<2) (|Gamma(VSWR=2)| = -9.54 dB) from computed Gamma
[~, fIdx] = min(abs(20*log10(abs(Gamma))+9.54));
BW_numeric = abs(f_r-f_meas(fIdx))/(f_r+f_meas(fIdx))*100;
fprintf('Numerically:  BW(VSWR<2) = %.2f %%\n', BW_numeric)

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
    B = W./(120*lambda0).*(1-0.636*log(k0*h));          % for h/lambda0 < 1/10
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
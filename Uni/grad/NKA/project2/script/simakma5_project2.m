%% Initialization
close all; clear; clc; addpath(genpath(fullfile(pwd, 'Uni', 'grad', 'NKA', 'project2', 'script')));
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";
fr1 = 0.9;
fr2 = 1.8;
Z0 = 50;

%% PIFA measurement reflection coefficient (linear)
PIFA_meas = sparameters(fullfile('data', 'pifa_sparam_measurement.S1P'));
f_PIFA_meas = PIFA_meas.Frequencies;
Gamma_PIFA_meas = squeeze(PIFA_meas.Parameters(1,1,:));

figure('Name', 'PIFA measurement reflection coefficient (linear)', 'Units', 'normalized', 'Position', [0, 0, 1, 1])
rfplot(PIFA_meas);
hold on
line = xline(fr1, '--', num2str(fr1));
line.LabelVerticalAlignment = 'top';
line = xline(fr2, '--', num2str(fr2));
line.LabelVerticalAlignment = 'top';
hold off
xlim([f_PIFA_meas(1)*1e-9, f_PIFA_meas(end)*1e-9])
ylim([-20, 0])
grid on
grid minor
xlabel('Frequency [GHz]')
ylabel('|\Gamma| [dB]')
legend('hide')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 26)
saveas(gcf, fullfile(pwd, '\Uni', '\grad', '\NKA', '\project2', '\report', '\src', '\pifa-meas-reflection-linear'), 'epsc')

%% PIFA reflection coefficient (linear)
PIFA_sim = sparameters(fullfile('data', 'pifa_sparam_simulation.S1P'));
f_PIFA_sim = PIFA_sim.Frequencies;
Gamma_PIFA_sim = squeeze(PIFA_sim.Parameters(1,1,:));

figure('Name', 'PIFA reflection coefficient (linear)', 'Units', 'normalized', 'Position', [0, 0, 1, 1])
rfplot(PIFA_sim);
hold on
rfplot(PIFA_meas);
line = xline(fr1, '--', num2str(fr1));
line.LabelVerticalAlignment = 'top';
line = xline(fr2, '--', num2str(fr2));
line.LabelVerticalAlignment = 'top';
line = yline(-12, '--', num2str(-12));
line.LabelVerticalAlignment = 'top';
hold off
xlim([max(f_PIFA_meas(1), f_PIFA_sim(1))*1e-9, min(f_PIFA_meas(end), f_PIFA_sim(end))*1e-9])
ylim([-20, 0])
grid on
grid minor
xlabel('Frequency [GHz]')
ylabel('|\Gamma| [dB]')
legend('Simulation', 'Measurement', 'Location', 'best')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 26)
saveas(gcf, fullfile(pwd, '\Uni', '\grad', '\NKA', '\project2', '\report', '\src', '\pifa-reflection-linear'), 'epsc')

%% PIFA reflection coefficient (Smith chart)
figure('Name', 'PIFA reflection coefficient (Smith chart)');
smithplot(PIFA_sim, 'GridBackgroundColor', 'w')
hold on
smith = smithplot(PIFA_meas, 'GridBackgroundColor', 'w');
hold off
smith.LegendLabels = {'Simulation', 'Measurement'};
saveas(gcf, fullfile(pwd, '\Uni', '\grad', '\NKA', '\project2', '\report', '\src', '\pifa-reflection-smith'), 'epsc')

%% PIFA input impedance
Zin_PIFA_sim = 50*(1+Gamma_PIFA_sim)./(1-Gamma_PIFA_sim);
Zin_PIFA_meas = 50*(1+Gamma_PIFA_meas)./(1-Gamma_PIFA_meas);

% figure('Name', 'PIFA input impedance', 'Units', 'normalized', 'Position', [0, 0, 1, 1])
plot(f_PIFA_sim*1e-9, real(Zin_PIFA_sim));
hold on
plot(f_PIFA_sim*1e-9, imag(Zin_PIFA_sim));
plot(f_PIFA_meas*1e-9, real(Zin_PIFA_meas));
plot(f_PIFA_meas*1e-9, imag(Zin_PIFA_meas));
line = xline(fr1, '--', num2str(fr1));
line.LabelVerticalAlignment = 'top';
line = xline(fr2, '--', num2str(fr2));
line.LabelVerticalAlignment = 'top';
line = yline(Z0, '--', num2str(Z0));
line.LabelVerticalAlignment = 'top';
line = yline(0, '--', num2str(0));
line.LabelVerticalAlignment = 'top';
hold off
xlim([max(f_PIFA_meas(1), f_PIFA_sim(1))*1e-9, min(f_PIFA_meas(end), f_PIFA_sim(end))*1e-9])
grid on
grid minor
xlabel('Frequency [GHz]')
ylabel('Z_{in} [\Omega]')
legend('R_{in} (sim.)', 'X_{in} (sim.)', 'R_{in} (meas.)', 'X_{in} (meas.)', 'location', 'best')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 26)
saveas(gcf, fullfile(pwd, '\Uni', '\grad', '\NKA', '\project2', '\report', '\src', '\pifa-impedance'), 'epsc')

%% PIFA radiation pattern (E-cut, 0.9 GHz)
% fid = fopen(fullfile('data', 'pifa_e-cut_0G9Hz.txt'), 'rt');
% Data = textscan(fid, '%f %f %f %f %f %f %f %f', 'headerLines', 4, 'CollectOutput', true);
% fclose(fid);
% sim = Data{1,1};
% theta_sim = sim(:,1)*pi/180;                        % convert to rad
% f_PIFA_sim = sim(sim(:,2)==90,3)-2.15;                   % convert from dBi to dB and
% f_PIFA_sim = [f_PIFA_sim; -(sim(sim(:,2)==270,3)-2.15)];      % flip data for theta==270

meas = readtable("pifa_e_cut_pattern.xlsx");
theta_meas = meas{3:end,1}*pi/180;                  % convert to rad
F_PIFA_meas = meas{3:end,3};

figure('Name', 'PIFA radiation pattern in the E-plane, 0.9 GHz');
polarplot(theta_meas, F_PIFA_meas)
rlim([-15, 0])
legend('hide')
saveas(gcf, fullfile(pwd, '\Uni', '\grad', '\NKA', '\project2', '\report', '\src', '\pifa-meas-radiation-e-0G9Hz'), 'epsc')

%% PIFA radiation pattern (H-cut, 0.9 GHz)
% fid = fopen(fullfile('data', 'pifa_h-cut_0G9Hz.txt'), 'rt');
% Data = textscan(fid, '%f %f %f %f %f %f %f %f', 'headerLines', 4, 'CollectOutput', true);
% fclose(fid);
% sim = Data{1,1};
% theta_sim = sim(:,1)*pi/180;                        % convert to rad
% f_PIFA_sim = sim(sim(:,2)==0,3)-2.15;                    % convert from dBi to dB and
% f_PIFA_sim = [f_PIFA_sim; -(sim(sim(:,2)==180,3)-2.15)];      % flip data for theta==270

meas = readtable("pifa_h_cut_pattern.xlsx");
theta_meas = meas{3:end,1}*pi/180;                  % convert to rad
F_PIFA_meas = meas{3:end,3};

figure('Name', 'PIFA radiation pattern in the H-plane, 0.9 GHz');
polarplot(theta_meas, F_PIFA_meas)
rlim([-15, 0])
legend('hide')
saveas(gcf, fullfile(pwd, '\Uni', '\grad', '\NKA', '\project2', '\report', '\src', '\pifa-meas-radiation-h-0G9Hz'), 'epsc')

%% PIFA radiation pattern (E-cut, 1.8 GHz)
meas = readtable("pifa_e_cut_pattern.xlsx");
theta_meas = meas{3:end,1}*pi/180;                  % convert to rad
F_PIFA_meas = meas{3:end,21};

figure('Name', 'PIFA radiation pattern in the E-plane, 0.9 GHz');
polarplot(theta_meas, F_PIFA_meas)
rlim([-25, 0])
legend('hide')
saveas(gcf, fullfile(pwd, '\Uni', '\grad', '\NKA', '\project2', '\report', '\src', '\pifa-meas-radiation-e-1G8Hz'), 'epsc')

%% PIFA radiation pattern (H-cut, 1.8 GHz)
meas = readtable("pifa_h_cut_pattern.xlsx");
theta_meas = meas{3:end,1}*pi/180;                  % convert to rad
F_PIFA_meas = meas{3:end,21};

figure('Name', 'PIFA radiation pattern in the H-plane, 0.9 GHz');
polarplot(theta_meas, F_PIFA_meas)
rlim([-25, 0])
legend('hide')
saveas(gcf, fullfile(pwd, '\Uni', '\grad', '\NKA', '\project2', '\report', '\src', '\pifa-meas-radiation-h-1G8Hz'), 'epsc')

%% Lambda-tenth monopole reflection coefficient (linear)
lambda10 = sparameters(fullfile('data', 'lambda_tenth_sparam_simulation.S1P'));
f_lambda10 = lambda10.Frequencies;
Gamma_lambda10 = squeeze(lambda10.Parameters(1,1,:));

figure('Name', 'Lambda-tenth reflection coefficient (linear)', 'Units', 'normalized', 'Position', [0, 0, 1, 1])
rfplot(lambda10);
line = xline(fr2, '--', num2str(fr2));
line.LabelVerticalAlignment = 'top';
line = yline(-15, '--', num2str(-15));
line.LabelVerticalAlignment = 'top';
hold off
xlim([1.3, 2.3])
ylim([-25, 0])
grid on
grid minor
xlabel('Frequency [GHz]')
ylabel('|\Gamma| [dB]')
legend('hide')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 26)
saveas(gcf, fullfile(pwd, '\Uni', '\grad', '\NKA', '\project2', '\report', '\src', '\lambda-tenth-reflection-linear'), 'epsc')

%% Lambda-tenth monopole reflection coefficient (Smith chart)
figure('Name', 'Lambda-tenth reflection coefficient (Smith chart)');
smithplot(lambda10, 'GridBackgroundColor', 'w')
legend('hide')
saveas(gcf, fullfile(pwd, '\Uni', '\grad', '\NKA', '\project2', '\report', '\src', '\lambda-tenth-reflection-smith'), 'epsc')

%% Lambda-tenth monopole input impedance
Zin_lambda10 = 50*(1+Gamma_lambda10)./(1-Gamma_lambda10);

figure('Name', 'Lambda-tenth input impedance', 'Units', 'normalized', 'Position', [0, 0, 1, 1])
plot(f_lambda10*1e-9, real(Zin_lambda10));
hold on
plot(f_lambda10*1e-9, imag(Zin_lambda10));
line = xline(fr2, '--', num2str(fr2));
line.LabelVerticalAlignment = 'top';
line = yline(Z0, '--', num2str(Z0));
line.LabelVerticalAlignment = 'top';
line = yline(0, '--', num2str(0));
line.LabelVerticalAlignment = 'top';
hold off
xlim([1.3, 2.3])
grid on
grid minor
xlabel('Frequency [GHz]')
ylabel('Z_{in} [\Omega]')
legend('R_{in} (sim.)', 'X_{in} (sim.)', 'location', 'best')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 26)
saveas(gcf, fullfile(pwd, '\Uni', '\grad', '\NKA', '\project2', '\report', '\src', '\lambda-tenth-impedance'), 'epsc')

%% Lambda-twentieth monopole reflection coefficient (linear)
lambda20 = sparameters(fullfile('data', 'lambda_twentieth_sparam_simulation.S1P'));
f_lambda20 = lambda20.Frequencies;
Gamma_lambda20 = squeeze(lambda20.Parameters(1,1,:));

figure('Name', 'Lambda-twentieth reflection coefficient (linear)', 'Units', 'normalized', 'Position', [0, 0, 1, 1])
rfplot(lambda20);
line = xline(fr2, '--', num2str(fr2));
line.LabelVerticalAlignment = 'top';
line = yline(-15, '--', num2str(-15));
line.LabelVerticalAlignment = 'top';
hold off
xlim([1.3, 2.3])
ylim([-20, 0])
grid on
grid minor
xlabel('Frequency [GHz]')
ylabel('|\Gamma| [dB]')
legend('hide')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 26)
saveas(gcf, fullfile(pwd, '\Uni', '\grad', '\NKA', '\project2', '\report', '\src', '\lambda-twentieth-reflection-linear'), 'epsc')

%% Lambda-twentieth monopole reflection coefficient (Smith chart)
figure('Name', 'Lambda-twentieth reflection coefficient (Smith chart)');
smithplot(lambda20, 'GridBackgroundColor', 'w')
legend('hide')
saveas(gcf, fullfile(pwd, '\Uni', '\grad', '\NKA', '\project2', '\report', '\src', '\lambda-twentieth-reflection-smith'), 'epsc')

%% Lambda-twentieth monopole input impedance
Zin_lambda20 = 50*(1+Gamma_lambda20)./(1-Gamma_lambda20);

figure('Name', 'Lambda-twentieth input impedance', 'Units', 'normalized', 'Position', [0, 0, 1, 1])
plot(f_lambda20*1e-9, real(Zin_lambda20));
hold on
plot(f_lambda20*1e-9, imag(Zin_lambda20));
line = xline(fr2, '--', num2str(fr2));
line.LabelVerticalAlignment = 'top';
line = yline(Z0, '--', num2str(Z0));
line.LabelVerticalAlignment = 'top';
line = yline(0, '--', num2str(0));
line.LabelVerticalAlignment = 'top';
hold off
xlim([1.3, 2.3])
grid on
grid minor
xlabel('Frequency [GHz]')
ylabel('Z_{in} [\Omega]')
legend('R_{in} (sim.)', 'X_{in} (sim.)', 'location', 'best')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 26)
saveas(gcf, fullfile(pwd, '\Uni', '\grad', '\NKA', '\project2', '\report', '\src', '\lambda-twentieth-impedance'), 'epsc')

%% Quality factor Q_3dB from Gamma drop of 3 dB
% PIFA
% - lower band
[Smin, frIdx] = min(20*log10(abs(PIFA_sim.Parameters(1,1,1:ceil(end/2)))));
fr = PIFA_sim.Frequencies(frIdx);
[~, f1Idx] = min(abs(20*log10(abs(PIFA_sim.Parameters(1,1,1:ceil(end/2))))-Smin-3));
f1 = PIFA_sim.Frequencies(f1Idx);
if f1 < fr
    f2 = PIFA_sim.Frequencies(f1Idx+2*(frIdx-f1Idx));
else
    f2 = f1;
    f1 = PIFA_sim.Frequencies(f1Idx-2*(f1Idx-frIdx));
end
Q_3dB_PIFA_LB = qdb_calc(f1, f2, fr, -3);

% - upper band
[Smin, frIdx] = min(20*log10(abs(PIFA_sim.Parameters(1,1,ceil(end/2):end))));
fr = PIFA_sim.Frequencies(ceil(end/2)+frIdx);
[~, f1Idx] = min(abs(20*log10(abs(PIFA_sim.Parameters(1,1,ceil(end/2):end)))-Smin-3));
f1 = PIFA_sim.Frequencies(ceil(end/2)+f1Idx);
if f1 < fr
    f2 = PIFA_sim.Frequencies(ceil(end/2)+f1Idx+2*(frIdx-f1Idx));
else
    f2 = f1;
    f1 = PIFA_sim.Frequencies(ceil(end/2)+f1Idx-2*(f1Idx-frIdx));
end
Q_3dB_PIFA_UB = qdb_calc(f1, f2, fr, -3);

% lambda-tenth monopole
[Smin, frIdx] = min(20*log10(abs(Gamma_lambda10(:))));
fr = f_lambda10(frIdx);
[~, f1Idx] = min(abs(20*log10(abs(Gamma_lambda10(:)))-Smin-3));
f1 = f_lambda10(f1Idx);
if f1 < fr
    f2 = f_lambda10(f1Idx+2*(frIdx-f1Idx));
else
    f2 = f1;
    f1 = f_lambda10(f1Idx-2*(f1Idx-frIdx));
end
Q_3dB_lambda10 = qdb_calc(f1, f2, fr, -3);

% lambda-twentieth monopole
[Smin, frIdx] = min(20*log10(abs(Gamma_lambda20(:))));
fr = f_lambda20(frIdx);
[~, f1Idx] = min(abs(20*log10(abs(Gamma_lambda20(:)))-Smin-3));
f1 = f_lambda20(f1Idx);
if f1 < fr
    f2 = f_lambda20(f1Idx+2*(frIdx-f1Idx));
else
    f2 = f1;
    f1 = f_lambda20(f1Idx-2*(f1Idx-frIdx));
end
Q_3dB_lambda20 = qdb_calc(f1, f2, fr, -3);

ka = 0.1:0.01:5;
Q_McLean = 1./(ka).^3 + 1./ka;
Q_Thal = 1.5./(ka).^3 + 0.6./ka;

ka_PIFA_LB = 2*pi*0.9e9/3e8*120e-3;     % largest dimension 120 mm (ground width)
ka_PIFA_UB = 2*pi*1.8e9/3e8*120e-3;     % largest dimension 120 mm (ground width)
ka_lambda10 = 2*pi*1.8e9/3e8*25e-3;     % largest dimension 25 mm (disc diameter)
ka_lambda20 = 2*pi*1.8e9/3e8*40e-3;     % largest dimension 40 mm (disc diameter)

[~, Q_Gustafsson, ka_Gustafsson, ~] = AntennaQ(linspace(0.0001,0.15,294), 0.006, 'rectangle_v', 2.45e9);
[~, kaIdx_Gustafsson_PIFA_LB] = min(abs(ka_Gustafsson-ka_PIFA_LB));
[~, kaIdx_Gustafsson_PIFA_UB] = min(abs(ka_Gustafsson-ka_PIFA_UB));
[~, kaIdx_Gustafsson_lambda10] = min(abs(ka_Gustafsson-ka_lambda10));
[~, kaIdx_Gustafsson_lambda20] = min(abs(ka_Gustafsson-ka_lambda20));

figure('Name', 'Quality factor from 3dB drop - fundamental limits', 'Units', 'normalized', 'Position', [0, 0, 1, 1]);
plot(ka, Q_McLean)
hold on
plot(ka, Q_Thal)
plot(ka_Gustafsson, Q_Gustafsson)
plot(ka_PIFA_LB, Q_3dB_PIFA_LB, 'x', 'MarkerSize', 16, 'LineWidth', 2)
plot(ka_PIFA_UB, Q_3dB_PIFA_UB, 'x', 'MarkerSize', 16, 'LineWidth', 2)
plot(ka_lambda10, Q_3dB_lambda10, 'x', 'MarkerSize', 16, 'LineWidth', 2)
plot(ka_lambda20, Q_3dB_lambda20, 'x', 'MarkerSize', 16, 'LineWidth', 2)
hold off
grid on
xlim([ka(1), ka(end)])
ylim([0, 1e3])
yscale(gca, 'log')
xlabel('ka [-]')
ylabel('Q [-]')
legend('Q_{McLean}', 'Q_{Thal}', 'Q_{Gustafsson}', 'Q_{3dB} - PIFA 0.9 GHz', ...
    'Q_{3dB} - PIFA 1.8 GHz', 'Q_{3dB} - lambda tenth', 'Q_{3dB} - lambda twentieth', ...
    'location', 'eastoutside')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 26)
saveas(gcf, fullfile(pwd, '\Uni', '\grad', '\NKA', '\project2', '\report', '\src', '\quality-factor-limits-q3db'), 'epsc')

%% Quality factor Q_Z from input impedance
% PIFA measurement
Rin = real(Zin_PIFA_meas);
Xin = imag(Zin_PIFA_meas);
dRin = derivace(2*pi*f_PIFA_meas, Rin);
dXin = derivace(2*pi*f_PIFA_meas, Xin);
dZin = abs(sqrt(dRin.^2 + (dXin + Xin./(2*pi*f_PIFA_meas)).^2));
Q_Z_PIFA_meas = (2*pi*f_PIFA_meas)./(2*Rin).*dZin;

% PIFA simulation
Rin = real(Zin_PIFA_sim);
Xin = imag(Zin_PIFA_sim);
dRin = derivace(2*pi*f_PIFA_sim, Rin);
dXin = derivace(2*pi*f_PIFA_sim, Xin);
dZin = abs(sqrt(dRin.^2 + (dXin + Xin./(2*pi*f_PIFA_sim)).^2));
Q_Z_PIFA_sim = (2*pi*f_PIFA_sim)./(2*Rin).*dZin;

% lambda-tenth monopole
Rin = real(Zin_lambda10);
Xin = imag(Zin_lambda10);
dRin = derivace(2*pi*f_lambda10, Rin);
dXin = derivace(2*pi*f_lambda10, Xin);
dZin = abs(sqrt(dRin.^2 + (dXin + Xin./(2*pi*f_lambda10)).^2));
Q_Z_lambda10 = (2*pi*f_lambda10)./(2*Rin).*dZin;

% lambda-twentieth monopole
Rin = real(Zin_lambda20);
Xin = imag(Zin_lambda20);
dRin = derivace(2*pi*f_lambda20, Rin);
dXin = derivace(2*pi*f_lambda20, Xin);
dZin = abs(sqrt(dRin.^2 + (dXin + Xin./(2*pi*f_lambda20)).^2));
Q_Z_lambda20 = (2*pi*f_lambda20)./(2*Rin).*dZin;

figure('Name', 'Quality factor from input impedance - fundamental limits', 'Units', 'normalized', 'Position', [0, 0, 1, 1]);
plot(ka, Q_McLean)
hold on
plot(ka, Q_Thal)
plot(ka_Gustafsson, Q_Gustafsson)
plot(ka_PIFA_LB, Q_Z_PIFA_sim(f_PIFA_sim == 0.9e9), 'x', 'MarkerSize', 16, 'LineWidth', 2)
plot(ka_PIFA_UB, Q_Z_PIFA_sim(f_PIFA_sim == 1.8e9), 'x', 'MarkerSize', 16, 'LineWidth', 2)
plot(ka_PIFA_LB, Q_Z_PIFA_meas(f_PIFA_meas == 0.9e9), 'x', 'MarkerSize', 16, 'LineWidth', 2)
plot(ka_PIFA_UB, Q_Z_PIFA_meas(f_PIFA_meas == 1.8e9), 'x', 'MarkerSize', 16, 'LineWidth', 2)
plot(ka_lambda10, Q_Z_lambda10(f_lambda10 == 1.8e9), 'x', 'MarkerSize', 16, 'LineWidth', 2)
plot(ka_lambda20, Q_Z_lambda20(f_lambda20 == 1.8e9), 'x', 'MarkerSize', 16, 'LineWidth', 2)
hold off
grid on
xlim([ka(1), ka(end)])
ylim([0, 1e3])
yscale(gca, 'log')
xlabel('ka [-]')
ylabel('Q [-]')
legend('Q_{McLean}', 'Q_{Thal}', 'Q_{Gustafsson}', 'Q_{Z} - PIFA sim. 0.9 GHz', ...
    'Q_{Z} - PIFA sim. 1.8 GHz', 'Q_{Z} - PIFA meas. 0.9 GHz', ...
    'Q_{Z} - PIFA meas. 1.8 GHz', 'Q_{Z} - lambda tenth', 'Q_{Z} - lambda twentieth', ...
    'location', 'eastoutside')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 26)
saveas(gcf, fullfile(pwd, '\Uni', '\grad', '\NKA', '\project2', '\report', '\src', '\quality-factor-limits-qz'), 'epsc')

%%
fprintf('\n\nSimulated PIFA lower band (0.9 GHz):\n----\n')
fprintf('ka:           %.2f\n', ka_PIFA_LB)
fprintf('Q_3dB:        %.2f\n', Q_3dB_PIFA_LB)
fprintf('Q_3dB/Q_min:  %.2f\n', Q_3dB_PIFA_LB/Q_Gustafsson(kaIdx_Gustafsson_PIFA_LB))
fprintf('Q_Z:          %.2f\n', Q_Z_PIFA_sim(f_PIFA_sim == 0.9e9))
fprintf('Q_Z/Q_min:    %.2f\n', Q_Z_PIFA_sim(f_PIFA_sim == 0.9e9)/Q_Gustafsson(kaIdx_Gustafsson_PIFA_LB))
fprintf('BW_3dB:       %.2f\n', bw_calc(Q_3dB_PIFA_LB,-3))
fprintf('BW_Q_Z:       %2.f\n', bw_calc(Q_Z_PIFA_sim(f_PIFA_sim == 0.9e9),-3))

fprintf('\nSimulated PIFA upper band (1.8 GHz):\n----\n')
fprintf('ka:           %.2f\n', ka_PIFA_UB)
fprintf('Q_3dB:        %.2f\n', Q_3dB_PIFA_UB)
fprintf('Q_3dB/Q_min:  %.2f\n', Q_3dB_PIFA_UB/Q_Gustafsson(kaIdx_Gustafsson_PIFA_UB))
fprintf('Q_Z:          %.2f\n', Q_Z_PIFA_sim(f_PIFA_sim == 1.8e9))
fprintf('Q_Z/Q_min:    %.2f\n', Q_Z_PIFA_sim(f_PIFA_sim == 1.8e9)/Q_Gustafsson(kaIdx_Gustafsson_PIFA_UB))
fprintf('BW_3dB:       %.2f\n', bw_calc(Q_3dB_PIFA_UB,-3))
fprintf('BW_Q_Z:       %2.f\n', bw_calc(Q_Z_PIFA_sim(f_PIFA_sim == 1.8e9),-3))

fprintf('\nMeasured PIFA lower band (0.9 GHz):\n----\n')
fprintf('ka:           %.2f\n', ka_PIFA_LB)
fprintf('Q_3dB:        ---\n')
fprintf('Q_3dB/Q_min:  ---\n')
fprintf('Q_Z:          %.2f\n', Q_Z_PIFA_meas(f_PIFA_sim == 0.9e9))
fprintf('Q_Z/Q_min:    %.2f\n', Q_Z_PIFA_meas(f_PIFA_sim == 0.9e9)/Q_Gustafsson(kaIdx_Gustafsson_PIFA_LB))
fprintf('BW_3dB:       ---\n')
fprintf('BW_Q_Z:       %2.f\n', bw_calc(Q_Z_PIFA_meas(f_PIFA_sim == 0.9e9),-3))

fprintf('\nMeasured PIFA upper band (1.8 GHz):\n----\n')
fprintf('ka:           %.2f\n', ka_PIFA_UB)
fprintf('Q_3dB:        ---\n')
fprintf('Q_3dB/Q_min:  ---\n')
fprintf('Q_Z:          %.2f\n', Q_Z_PIFA_meas(f_PIFA_sim == 1.8e9))
fprintf('Q_Z/Q_min:    %.2f\n', Q_Z_PIFA_meas(f_PIFA_sim == 1.8e9)/Q_Gustafsson(kaIdx_Gustafsson_PIFA_UB))
fprintf('BW_3dB:       ---\n')
fprintf('BW_Q_Z:       %2.f\n', bw_calc(Q_Z_PIFA_meas(f_PIFA_sim == 1.8e9),-3))

fprintf('\nLambda-tenth monopole:\n----\n')
fprintf('ka:           %.2f\n', ka_lambda10)
fprintf('Q_3dB:        %.2f\n', Q_3dB_lambda10)
fprintf('Q_3dB/Q_min:  %.2f\n', Q_3dB_lambda10/Q_Gustafsson(kaIdx_Gustafsson_lambda10))
fprintf('Q_Z:          %.2f\n', Q_Z_lambda10(f_lambda10 == 1.8e9))
fprintf('Q_Z/Q_min:    %.2f\n', Q_Z_lambda10(f_lambda10 == 1.8e9)/Q_Gustafsson(kaIdx_Gustafsson_lambda10))
fprintf('BW_3dB:       %.2f\n', bw_calc(Q_3dB_lambda10,-3))
fprintf('BW_Q_Z:       %2.f\n', bw_calc(Q_Z_lambda10(f_lambda10 == 1.8e9),-3))

fprintf('\nLambda-twentieth monopole:\n----\n')
fprintf('ka:           %.2f\n', ka_lambda20)
fprintf('Q_3dB:        %.2f\n', Q_3dB_lambda20)
fprintf('Q_3dB/Q_min:  %.2f\n', Q_3dB_lambda20/Q_Gustafsson(kaIdx_Gustafsson_lambda20))
fprintf('Q_Z:          %.2f\n', Q_Z_lambda20(f_lambda20 == 1.8e9))
fprintf('Q_Z/Q_min:    %.2f\n', Q_Z_lambda20(f_lambda20 == 1.8e9)/Q_Gustafsson(kaIdx_Gustafsson_lambda20))
fprintf('BW_3dB:       %.2f\n', bw_calc(Q_3dB_lambda20,-3))
fprintf('BW_Q_Z:       %2.f\n', bw_calc(Q_Z_lambda20(f_lambda20 == 1.8e9),-3))

%% Functions
function qdb = qdb_calc(f1, f2, fc, db)
    gam = db2mag(db);
    psv = (1+gam)./(1-gam);
    bw = (f2 - f1)./fc;

    qdb = (psv - 1)./(sqrt(psv)*bw)/10;
end

function bw = bw_calc(q, db) % bw jednotka -> %
    if nargin < 2
        db = -3;
    end

    gam = db2mag(db);
    psv = (1+gam)./(1-gam);

    bw = (psv - 1)./(q*sqrt(psv))*100;
end

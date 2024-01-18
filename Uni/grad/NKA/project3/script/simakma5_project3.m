%% Initialization
close all; clear; clc; addpath(genpath(fullfile(pwd, 'Uni', 'grad', 'NKA', 'project3', 'script')));
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";
fr = 2.45;
Z0 = 50;

%% Reflection coefficient
sim = sparameters(fullfile('data', 'reflection-simulation.S1P'));
f_sim = sim.Frequencies;
Gamma_sim = squeeze(sim.Parameters(1,1,:));

meas = sparameters(fullfile('data', 'reflection-measurement.S1P'));
f_meas = meas.Frequencies;
Gamma_meas = squeeze(meas.Parameters(1,1,:));

figure('Name', 'Reflection coefficint', 'Units', 'normalized', 'Position', [0, 0, 1, 1])
rfplot(sim)
hold on
rfplot(meas)
xline(fr, '--', num2str(fr))
yline(-10, '--', num2str(-10))
hold off
grid on
grid minor
xlim([max(f_sim(1), f_meas(1))*1e-9, min(f_sim(end), f_meas(end))*1e-9])
xlabel('Frequency [GHz]')
ylabel('|\Gamma| [dB]')
legend('Simulation', 'Measurement', 'location', 'southwest')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 26)
saveas(gcf, fullfile(pwd, '\Uni', '\grad', '\NKA', '\project3', '\report', '\src', '\reflection-coefficient'), 'epsc')

%% Input impedance
Zin_sim = 50*(1+Gamma_sim)./(1-Gamma_sim);
Zin_meas = 50*(1+Gamma_meas)./(1-Gamma_meas);

figure('Name', 'Input impedance', 'Units', 'normalized', 'Position', [0, 0, 1, 1])
plot(f_sim*1e-9, real(Zin_sim))
hold on
plot(f_sim*1e-9, imag(Zin_sim))
plot(f_meas*1e-9, real(Zin_meas))
plot(f_meas*1e-9, imag(Zin_meas))
xline(fr, '--', num2str(fr))
yline(50, '--', num2str(50))
hold off
grid on
grid minor
xlim([max(f_sim(1), f_meas(1))*1e-9, min(f_sim(end), f_meas(end))*1e-9])
xlabel('Frequency [GHz]')
ylabel('Input impedance [\Omega]')
legend('R_{in} (sim.)', 'X_{in} (sim.)', 'R_{in} (meas.)', 'X_{in} (meas.)', 'location', 'northwest')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 26)
saveas(gcf, fullfile(pwd, '\Uni', '\grad', '\NKA', '\project3', '\report', '\src', '\input-impedance'), 'epsc')

%% Axial ratio
fid = fopen(fullfile('data', 'axial-ratio.txt'), 'rt');
Data = textscan(fid, '%f %f', 'headerLines', 3, 'CollectOutput', true);
fclose(fid);
AR = Data{1,1};
f_AR = AR(:,1);
AR = AR(:,2);

figure('Name', 'Axial ratio', 'Units', 'normalized', 'Position', [0, 0, 1, 1])
plot(f_AR, AR)
hold on
xline(fr, '--', num2str(fr))
axis tight
grid on
grid minor
ylabel('Axial ratio [dB]')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 26)
saveas(gcf, fullfile(pwd, '\Uni', '\grad', '\NKA', '\project3', '\report', '\src', '\axial-ratio'), 'epsc')

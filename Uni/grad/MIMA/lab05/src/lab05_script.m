%%
close all; clear; clc
% Read ENR values of the used noise source
T1 = table2array(readtable('ENR.txt', 'NumHeaderLines', 7));
f_ENR = T1(:,1)/1e6;
ENR = 10.^(T1(:,2)/10);

%% Task 1
% Read the measured noise power data
T1 = table2array(readtable('N1_OFF.TXT', 'NumHeaderLines', 28));
f = T1(:, 1)/1e6;
N1_COLD = 10.^(T1(:, 2)/10)*1e-3;
T1 = table2array(readtable('N1_ON.TXT', 'NumHeaderLines', 28));
% f = T1(:, 1)/1e6;
N1_HOT = 10.^(T1(:, 2)/10)*1e-3;
T1 = table2array(readtable('N2_OFF.TXT', 'NumHeaderLines', 28));
% f = T1(:, 1)/1e6;
N2_COLD = 10.^(T1(:, 2)/10)*1e-3;
T1 = table2array(readtable('N2_ON.TXT', 'NumHeaderLines', 28));
% f = T1(:, 1)/1e6;
N2_HOT = 10.^(T1(:, 2)/10)*1e-3;
% Plot
fig = figure(1);
plot( ...
    f,10*log10(N1_COLD/1e-3), ...
    f,10*log10(N1_HOT/1e-3), ...
    f,10*log10(N2_COLD/1e-3), ...
    f,10*log10(N2_HOT/1e-3) ...
    )
xlim([20 1.5e3])
xlabel('Frequency [MHz]')
ylabel('Noise power [dBm]')
grid on
grid minor
legend('N_{1,COLD}','N_{1,HOT}','N_{2,COLD}','N_{2,HOT}','Location','east')
% draw a white rectangle around the graph to avoid trimming during export
a = annotation("rectangle",[0 0 1 1],"Color",'w');
exportgraphics(fig,"task1_powers.eps")
delete(a)

% Calculate the noise figures
F1 = zeros(length(N1_COLD),1);
F2 = zeros(length(N1_COLD),1);
for i=1:length(f)
    [~,k] = min(abs(f_ENR-f(i)));
    F1(i) = ENR(k)/(N1_HOT(i)/N1_COLD(i)+1);
    F2(i) = ENR(k)/(N2_HOT(i)/N2_COLD(i)+1);
end
% Plot results
fig = figure(2);
plot( ...
    f,10*log10(F1), ...
    f,10*log10(F2) ...
    )
xlim([20 1.5e3])
xlabel('Frequency [MHz]')
ylabel('Noise figure [dB]')
grid on
grid minor
legend('F_1','F_2','Location','east')
% draw a white rectangle around the graph to avoid trimming during export
a = annotation("rectangle",[0 0 1 1],"Color",'w');
exportgraphics(fig,"task1-figures.eps")
delete(a)

% Plot equivalent noise temperatures
fig = figure(3);
plot( ...
    f,290*(F1-1), ...
    f,290*(F2-1) ...
    )
xlim([20 1.5e3])
xlabel('Frequency [MHz]')
ylabel('Equivalent noise temperature [K]')
grid on
grid minor
legend('T_1','T_2','Location','east')
% draw a white rectangle around the graph to avoid trimming during export
a = annotation("rectangle",[0 0 1 1],"Color",'w');
exportgraphics(fig,"task1-temperatures.eps")
delete(a)

%% Task 2
% Read the measured noise power data
T1 = table2array(readtable('N2_Y-AMP_OFF.TXT', 'NumHeaderLines', 28));
f = T1(:, 1)/1e6;
N2_COLD_AMP = 10.^(T1(:, 2)/10)*1e-3;
T1 = table2array(readtable('N2_Y-AMP_ON.TXT', 'NumHeaderLines', 28));
% f = T1(:, 1)/1e6;
N2_HOT_AMP = 10.^(T1(:, 2)/10)*1e-3;
T1 = table2array(readtable('N2_Y-ATT_OFF.TXT', 'NumHeaderLines', 28));
% f = T1(:, 1)/1e6;
N2_COLD_ATT = 10.^(T1(:, 2)/10)*1e-3;
T1 = table2array(readtable('N2_Y-ATT_ON.TXT', 'NumHeaderLines', 28));
% f = T1(:, 1)/1e6;
N2_HOT_ATT = 10.^(T1(:, 2)/10)*1e-3;
% Plot
fig = figure(4);
plot( ...
    f,10*log10(N2_COLD_AMP/1e-3), ...
    f,10*log10(N2_HOT_AMP/1e-3), ...
    f,10*log10(N2_COLD_ATT/1e-3), ...
    f,10*log10(N2_HOT_ATT/1e-3) ...
    )
xlim([20 1.5e3])
xlabel('Frequency [MHz]')
ylabel('Noise power [dBm]')
grid on
grid minor
legend({'N_{2,COLD}^{AMP}','N_{2,HOT}^{AMP}','N_{2,COLD}^{ATT}','N_{2,HOT}^{ATT}'},'Location','north')
% draw a white rectangle around the graph to avoid trimming during export
a = annotation("rectangle",[0 0 1 1],"Color",'w');
exportgraphics(fig,"task2-powers.eps")
delete(a)

% Calculate the gain of the DUTs
G_AMP = (N2_HOT_AMP-N2_COLD_AMP)./(N2_HOT-N2_COLD);
G_ATT = (N2_HOT_ATT-N2_COLD_ATT)./(N2_HOT-N2_COLD);
% Plot results
fig = figure(5);
plot( ...
    f,10*real(log10(G_AMP)), ...
    f,10*real(log10(G_ATT)) ...
    )
xlim([20 1.5e3])
xlabel('Frequency [MHz]')
ylabel('Gain [dB]')
grid on
grid minor
legend({'G^{AMP}','G^{ATT}'},'Location','Southeast')
% draw a white rectangle around the graph to avoid trimming during export
a = annotation("rectangle",[0 0 1 1],"Color",'w');
exportgraphics(fig,"task2-gains.eps")
delete(a)

% Calculate the noise figures
% First overall SPA + DUT
F_AMP = zeros(length(N2_COLD_AMP),1);
F_ATT = zeros(length(N2_COLD_AMP),1);
for i=1:length(f)
    [~,k] = min(abs(f_ENR-f(i)));
    F_AMP(i) = ENR(k)/(N2_HOT_AMP(i)/N2_COLD_AMP(i)-1);
    F_ATT(i) = ENR(k)/(N2_HOT_ATT(i)/N2_COLD_ATT(i)-1);
end
% Separate DUT from SPA using Frii's formula
F_AMP = F_AMP-(F2-1)./G_AMP;
F_ATT = F_ATT-(F2-1)./G_ATT;
% Plot results
fig = figure(6);
plot( ...
    f,10*real(log10(F_AMP)), ...
    f,10*real(log10(F_ATT)) ...
    )
xlim([20 1.5e3])
xlabel('Frequency [MHz]')
ylabel('Noise figure [dB]')
grid on
grid minor
legend('F^{AMP}','F^{ATT}','Location','Northeast')
% draw a white rectangle around the graph to avoid trimming during export
a = annotation("rectangle",[0 0 1 1],"Color",'w');
exportgraphics(fig,"task2-figures.eps")
delete(a)

% Plot equivalent noise temperatures
fig = figure(7);
plot( ...
    f,290*(F_AMP-1), ...
    f,290*(F_ATT-1) ...
    )
xlim([20 1.5e3])
ylim([-1e4 5e4])
xlabel('Frequency [MHz]')
ylabel('Equivalent noise temperature [K]')
grid on
grid minor
legend('T^{AMP}','T^{ATT}','Location','northeast')
% draw a white rectangle around the graph to avoid trimming during export
a = annotation("rectangle",[0 0 1 1],"Color",'w');
exportgraphics(fig,"task2-temperatures.eps")
delete(a)

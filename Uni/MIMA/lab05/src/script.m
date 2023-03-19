%%
close all; clear; clc

%% Task 1
clear
T1 = table2array(readtable('N1_OFF.TXT', 'NumHeaderLines', 28));
f = T1(:, 1)/1e6;
N = 10.^(T1(:, 2)/10)*1e-3;

figure(1)
plot(f,N)
xlim([20 1.5e3])
grid on
grid minor
hold on

T1 = table2array(readtable('N1_ON.TXT', 'NumHeaderLines', 28));
f = T1(:, 1)/1e6;
N = 10.^(T1(:, 2)/10)*1e-3;
plot(f,N)

T1 = table2array(readtable('N2_OFF.TXT', 'NumHeaderLines', 28));
f = T1(:, 1)/1e6;
N = 10.^(T1(:, 2)/10)*1e-3;
plot(f,N)

T1 = table2array(readtable('N2_ON.TXT', 'NumHeaderLines', 28));
f = T1(:, 1)/1e6;
N = 10.^(T1(:, 2)/10)*1e-3;
plot(f,N)

legend('N1_{OFF}','N1_{ON}','N2_{OFF}','N2_{ON}')
hold off

%% Task 2
clear
T1 = table2array(readtable('N2_Y-AMP_OFF.TXT', 'NumHeaderLines', 28));
f = T1(:, 1)/1e6;
N = 10.^(T1(:, 2)/10)*1e-3;

figure(1)
plot(f,N)
xlim([20 1.5e3])
grid on
grid minor
hold on

T1 = table2array(readtable('N2_Y-AMP_ON.TXT', 'NumHeaderLines', 28));
f = T1(:, 1)/1e6;
N = 10.^(T1(:, 2)/10)*1e-3;
plot(f,N)

T1 = table2array(readtable('N2_Y-ATT_OFF.TXT', 'NumHeaderLines', 28));
f = T1(:, 1)/1e6;
N = 10.^(T1(:, 2)/10)*1e-3;
plot(f,N)

T1 = table2array(readtable('N2_Y-ATT_ON.TXT', 'NumHeaderLines', 28));
f = T1(:, 1)/1e6;
N = 10.^(T1(:, 2)/10)*1e-3;
plot(f,N)

legend({'N2_{OFF}^{AMP}','N2_{ON}^{AMP}','N2_{OFF}^{ATT}','N2_{ON}^{ATT}'},'Location','Northwest')
hold off

%% Task 3

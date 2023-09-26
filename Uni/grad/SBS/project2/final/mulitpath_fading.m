close all; clear; clc

%% import the data
data13 = importdata('RecTrace_013.csv');
data14 = importdata('RecTrace_014.csv');
data15 = importdata('RecTrace_015.csv');

%% Preprocessing
for i = 1:156
    dataline13(i) = max(data13(i,:));
end
for i = 1:1509
    dataline14(i) = max(data14(i,:));
end
for i = 1:337
    dataline15(i) = max(data15(i,:));
end

figure(9)
subplot(311)
plot(dataline13)
xlabel('Samples [-]')
ylabel('Received power [dBm]')
title('Nezastinene mereni')
subplot(312)
plot(dataline14)
xlabel('Samples [-]')
ylabel('Received power [dBm]')
title('Dynamicky zastinene mereni')
subplot(313)
plot(dataline15)
xlabel('Samples [-]')
ylabel('Received power [dBm]')
title('Uplne zastinene mereni')

%% Calculate empirical CDF
[f13,x13] = ecdf(dataline13(83:141));   % relevant data selected
x_norm13 = x13-median(x13);
[f14,x14] = ecdf(dataline14);
x_norm14 = x14-median(x14);
[f15,x15] = ecdf(dataline15);
x_norm15 = x15-median(x15);


%% Generate rayleigh distribution CDF
xr = linspace(0.01,4,1000);
yr = cdf('Rayleigh',xr,1);
kr = find(yr>0.499 & yr<0.501);
xr_norm = 20*log10(xr)-20*log10((xr(1,kr)));

%% Results
figure(1)
plot(xr_norm,yr,'b--')
hold on
set(gca,'yscale','log')
plot(x_norm13,f13,'--')
plot(x_norm14,f14,'--')
plot(x_norm15,f15,'--')
xline(0,'--','Color',[0.4 0.4 0.4])
yline(90e-2,'--','Color',[0.4 0.4 0.4])
yline(50e-2,'--','Color',[0.4 0.4 0.4])
yline(10e-2,'--','Color',[0.4 0.4 0.4])
yline(5e-2,'--','Color',[0.4 0.4 0.4])
set(gca,'yscale','log')
xlabel('Relative loss L [dB], normalized to 0 dB at 50%','FontSize',16)
ylabel('CDF(L) [-]','FontSize',16)
legend('Rayleigh','Nezastinene mereni','Dynamicky zastinene mereni','Uplne zastinene mereni','Percentily 90, 50, 10, 5','Location','southeast')
axis square
ylim([.009 1])
grid on

disp(['LoS: min = ' num2str(min(dataline13(83:141))) ...
    ', max = ' num2str(max(dataline13(83:141))) ...
    ', median = ' num2str(median(dataline13(83:141)))])
disp(['Dynamic shadowing: min = ' num2str(min(dataline14)) ...
    ', max = ' num2str(max(dataline13(83:141))) ...
    ', median = ' num2str(median(dataline13(83:141)))])
disp(['Static shadowing: min = ' num2str(min(dataline15)) ...
    ', max = ' num2str(max(dataline15)) ...
    ', median = ' num2str(median(dataline15))])



close all; clear; clc

%% ELSY SG2000
disp("=============")
disp("ELSY SG2000")
disp("-------------")
P_set = [10 6 2 -2 -6 -10 -14 -18 -22 -26];
P_C = [10.5 6.4 2.3 -1.7 -5.7 -9.7 -13.8 -17.8 -21.8 -26];
P_IM2 = [-21.2 -31.8 -42.4 -47.7 -51.7 -64 -68.4 -72 -76 -78.6];
P_IM3 = [-29.1 -43.5 -56.7 -65.6 -70 -79 -79 -78];

% find intercept points IP2 and IP3
start = 10;
stop = 50;
num = 1e3;
X = linspace(start,stop,num);
coeff_0 = polyfit(P_set,P_C,1);
line_0 = coeff_0(1)*X + coeff_0(2);

coeff_IM2 = polyfit(P_set(1:5),P_IM2(1:5),1);
line_IM2 = coeff_IM2(1)*X + coeff_IM2(2);
[~,IIP2] = min(abs(line_0-line_IM2));
OIP2 = line_IM2(IIP2);
disp(['IP2 = ' num2str(round(OIP2,1))])

coeff_IM3 = polyfit(P_set(1:5),P_IM3(1:5),1);
line_IM3 = coeff_IM3(1)*X + coeff_IM3(2);
[~,IIP3] = min(abs(line_0-line_IM3));
OIP3 = line_IM3(IIP3);
disp(['IP3 = ' num2str(round(OIP3,1))])

% plot
fig = figure(1);
plot(P_set,P_C,P_set,P_IM2)
hold on
plot(P_set(1:length(P_set)-2),P_IM3)
% plot(P_set,polyval(coeff_IM2,P_set),'--','Color',[0.5 0.5 0.5])
% plot(P_set,polyval(coeff_IM3,P_set),'--','Color',[0.5 0.5 0.5])
xlabel('P_{set} [dBm]')
ylabel('[dBm]')
axis('tight')
grid on
grid minor
legend(["P_C" "P_{IM2}" "P_{IM3}"],"Location","Northwest")
hold off

% export
% draw a white rectangle around the graph to avoid trimming during export
a = annotation("rectangle",[0 0 1 1],"Color",'w');
exportgraphics(fig,"ELSY-SG2000-intermodulation.eps")
delete(a)

disp("=============")

%% R&S SMF 100A
disp("=============")
disp("R&S SMF 100A")
disp("-------------")
P_set = [10 5 0 -5 -10 -15 -20 -25 -30];
P_C = [10.2 5.1 0.1 -4.8 -9.8 -14.8 -19.8 -24.8 -29.8];
P_IM2 = [-30.3 -40.7 -50.8 -51 -61.3 -60.9 -70.6 -70.9 -77.5];
P_IM3 = [-55.1 -69.8 -77.3];

% find intercept points IP2 and IP3
start = 10;
stop = 50;
num = 1e3;
X = linspace(start,stop,num);
coeff_0 = polyfit(P_set,P_C,1);
line_0 = coeff_0(1)*X + coeff_0(2);

coeff_IM2 = polyfit(P_set(1:3),P_IM2(1:3),1);
line_IM2 = coeff_IM2(1)*X + coeff_IM2(2);
[~,IIP2] = min(abs(line_0-line_IM2));
OIP2 = line_IM2(IIP2);
disp(['IP2 = ' num2str(round(OIP2,1))])

coeff_IM3 = polyfit(P_set(1:3),P_IM3,1);
line_IM3 = coeff_IM3(1)*X + coeff_IM3(2);
[~,IIP3] = min(abs(line_0-line_IM3));
OIP3 = line_IM3(IIP3);
disp(['IP3 = ' num2str(round(OIP3,1))])

% plot
fig = figure(2);
plot(P_set,P_C,P_set,P_IM2)
hold on
plot(P_set(1:3),P_IM3)
% plot(P_set,polyval(coeff_IM2,P_set),'--','Color',[0.5 0.5 0.5])
% plot(P_set,polyval(coeff_IM3,P_set),'--','Color',[0.5 0.5 0.5])
xlabel('P_{set} [dBm]')
ylabel('[dBm]')
axis('tight')
grid on
grid minor
legend(["P_C" "P_{IM2}" "P_{IM3}"],"Location","Northwest")
hold off

% export
% draw a white rectangle around the graph to avoid trimming during export
a = annotation("rectangle",[0 0 1 1],"Color",'w');
exportgraphics(fig,"RS-SMF-100A-intermodulation.eps")
delete(a)

disp("=============")

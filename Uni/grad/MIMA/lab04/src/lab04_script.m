close all; clear; clc

%% Task 1
P_LO = [-10 -8 -6 -4 -2 0 2 4 6 8 10];
P_IF = [-64 -60 -55 -49.3 -42.4 -35.3 -29.2 -24.9 -22.6 -21.5 -21];
P_IFLO = [-37.4 -35.4 -33.4 -31.6 -29.7 -28 -26.7 -26 -26.1 -27.1 -28.5];
L_c = round(-10 - P_IF,1)
O_IFLO = round(P_IF - P_IFLO,1)

%% Task 2
P_RF = [-10 -8 -6 -4 -2 0];
P_IF = [-20.5 -18.4 -16.4 -14.4 -12.4 -10.4];
P_3 = [-85 -80 -75 -70 -66 -64];
O_IF3 = round(P_IF - P_3,1)

% approximate linear behvaiour at the beginning
coeff_IF = polyfit(P_RF(1:5),P_IF(1:5),1);
line_IF = polyval(coeff_IF,P_RF);
coeff_3 = polyfit(P_RF(1:5),P_3(1:5),1);
line_3 = polyval(coeff_3,P_RF);

% find IP3
X_alt = linspace(-10,40,1e4);
line_IF_alt = coeff_IF(1)*X_alt + coeff_IF(2);
line_3_alt = coeff_3(1)*X_alt + coeff_3(2);
diff = line_IF_alt - line_3_alt;
[~, i] = min(abs(diff));
IP3_in = round(X_alt(i),1);
IP3_out = round(line_IF_alt(i),1);
% disp(['IP3 = [' num2str(IP3_in) ',' num2str(IP3_out) ']'])
disp(['IP3 = ' num2str(IP3_out)])

% plot
fig = figure(1);
plot(P_RF,P_IF,P_RF,P_3)
hold on
X_alt = linspace(-10,20);
plot(X_alt,coeff_IF(1)*X_alt+coeff_IF(2),'--',X_alt,coeff_3(1)*X_alt+coeff_3(2),'--','Color',[0.5 0.5 0.5])
grid on
grid minor
legend(["P_{IF}" "P_3"],"Location","Southeast")
xlabel("P_{in} [dBm]");
ylabel("P_{out} [dBm]");
hold off
% draw a white rectangle around the graph to avoid trimming during export
a = annotation("rectangle",[0 0 1 1],"Color",'w');
exportgraphics(fig,"task2-ip3.eps")
delete(a)

%% Task 3a
P_IF1 = [-21.1 -20.9 -21.2 -21.3];
P_IF2 = [-21.3 -21.6 -22.2 -23.9];
L_L = round(-10-P_IF1,1)
L_U = round(-10-P_IF2,1)
L_M = round(P_IF1-P_IF2,1)

%% Task 3b
P_IF1 = [-21.6 -20.8 -21.3 -21.1];
P_IF2 = [-40.9 -48 -42 -56.6];
L_L = round(-10-P_IF1,1)
L_U = round(-10-P_IF2,1)
L_M = round(P_IF1-P_IF2,1)

%% Task 4a
P_IN = [0 2 4 6 8 10];
P_OUT2 = [-49.3 -45.2 -41.2 -38.5 -36.8 -35.3];
P_OUT3 = [-49.2 -47.5 -45.8 -43.7 -43 -44.2];
P_OUT4 = [-55.1 -53.7 -52.9 -49.3 -43.2 -39.1];
L_k = round(P_IN-P_OUT2,1)
O_OUT = zeros(1,length(P_OUT2));
for i = 1:length(P_OUT2)
    O_OUT(i) = round(P_OUT2(i)-max(P_OUT3(i),P_OUT4(i)),1);
end
O_OUT

%% Task 4b
P_IN = [0 2 4 6 8 10];
P_OUT2 = [-88 -85 -85 -88 -79 -71.8];
P_OUT3 = [-45.6 -38.7 -30.8 -18.6 -9.8 -5.1];
P_OUT4 = [-88 -88 -88 -79 -79 -76];
L_k = round(P_IN-P_OUT2,1)
O_OUT = zeros(1,length(P_OUT2));
for i = 1:length(P_OUT2)
    O_OUT(i) = round(P_OUT2(i)-max(P_OUT3(i),P_OUT4(i)),1);
end
O_OUT

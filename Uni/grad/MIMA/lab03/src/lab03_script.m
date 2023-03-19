%%
close all; clear; clc;

%% Task 1: 1dB compression point measurement
% Amplifier ERA-3SM
L_att = [30 25 22 20 19 18 17 16 15 11 9 8 7 6 5];
P_out = [-10.5 -5.1 -2.3 -0.2 0.9 1.9 2.7 3.7 4.8 8.6 10.4 11.1 11.7 12.3 12.7];
% G = L_att + P_out; disp(G)

% approximate linear behvaiour at the beginning
coeff = polyfit(-L_att(1:8),P_out(1:8),1);
line = polyval(coeff,-L_att);

% find the point of 1dB compression
diff = line - P_out - 1;
[~, i] = min(abs(diff));
P1dB_in = -L_att(i);
P1dB_out = P_out(i);
disp(['ERA-3SM+: P_1dB = ' num2str(P1dB_out)])

% plot
fig = figure(1);
p = plot(-L_att,P_out,'MarkerIndices',i);
% xticks([-30 -25 -20 -15 -10 P1dB_in -5])
% yticks([-8 -4 0 4 8 P1dB_out])
% p.Marker = '*';
axis('tight')
hold on
p = plot(-L_att,line,'--','Color',[0.5 0.5 0.5],'MarkerIndices',i);
% p.Marker = '*';
xline(P1dB_in,'--','Color',[0.5 0.5 0.5])
% yline(P1dB_out,'--','Color',[0.5 0.5 0.5])
grid on
grid minor
% legend("boxoff")
xlabel("P_{in} [dBm]");
ylabel("P_{out} [dBm]");
% title("ERA-3SM")
hold off

% draw a white rectangle around the graph to avoid trimming during export
a = annotation("rectangle",[0 0 1 1],"Color",'w');
exportgraphics(fig,"ERA3SM-compression.eps")
delete(a)


% Amplifier BFU-760F
L_att = [30 25 20 18 15 13 12 11 10 9 8];
P_out = [-9.9 -4.5 0.4 2.4 5.3 7 7.8 8.8 9.5 10.1 10.5];
% G = L_att + P_out; disp(G)

% approximate linear behvaiour at the beginning
coeff = polyfit(-L_att(1:6),P_out(1:6),1);
line = polyval(coeff,-L_att);

% find the point of 1dB compression
diff = line - P_out - 1;
[~, i] = min(abs(diff));
P1dB_in = -L_att(i);
P1dB_out = P_out(i);
disp(['LNA with BFU-760F: P_1dB = ' num2str(P1dB_out)])

% plot
fig = figure(2);
p = plot(-L_att,P_out,'MarkerIndices',i);
% xticks([-30 -25 -20 -15 P1dB_in])
% yticks([-8 -4 0 4 8 P1dB_out])
% p.Marker = '*';
axis('tight')
hold on
p = plot(-L_att,line,'--','Color',[0.5 0.5 0.5],'MarkerIndices',i);
% p.Marker = '*';
xline(P1dB_in,'--','Color',[0.5 0.5 0.5])
% yline(P1dB_out,'--','Color',[0.5 0.5 0.5])
grid on
grid minor
% legend("boxoff")
xlabel("P_{in} [dBm]");
ylabel("P_{out} [dBm]");
% title("BFU-760F")
hold off

% draw a white rectangle around the graph to avoid trimming during export
a = annotation("rectangle",[0 0 1 1],"Color",'w');
exportgraphics(fig,"BFU760F-compression.eps")
delete(a)


%% Task 2: IP3 measurement
% Amplifier ERA-3SM
L_att = [23 20 17 15 14 11 8 7 5 2];
P_a = [-3.2 -0.2 2.7 4.6 5.5 7.6 8.9 9.2 9.6 10.1];
P_IM3 = [-60 -50.6 -40.7 -33.4 -29.7 -17.8 -10.6 -9.4 -7.5 -7.6];
% O_IM3 = P_a - P_IM3; disp(O_IM3)

% approximate linear behvaiour at the beginning
coeff_a = polyfit(-L_att(1:5),P_a(1:5),1);
line_a = polyval(coeff_a,-L_att);
coeff_IM3 = polyfit(-L_att(1:5),P_IM3(1:5),1);
line_IM3 = polyval(coeff_IM3,-L_att);

% find IP3
X_alt = linspace(-5,10);
line_a_alt = coeff_a(1)*X_alt + coeff_a(2);
line_IM3_alt = coeff_IM3(1)*X_alt + coeff_IM3(2);
diff = line_a_alt - line_IM3_alt;
[~, i] = min(abs(diff));
IP3_in = round(X_alt(i),1);
IP3_out = round(line_a_alt(i),1);
disp(['ERA-3SM+: IP3 = ' num2str(IP3_out)])

% plot
fig = figure(3);
plot(-L_att,P_a,-L_att,P_IM3)
axis('tight')
hold on
plot(-L_att,line_a,'--',-L_att,line_IM3,'--','Color',[0.5 0.5 0.5])
grid on
grid minor
legend(["P_a" "P_{2a \pm b}"],"Location","Southeast")
% legend("boxoff")
ax = gca;
xlabel("P_{in} [dBm]");
ylabel("P_{out} [dBm]");
% title("ERA-3SM")
hold off

% draw a white rectangle around the graph to avoid trimming during export
a = annotation("rectangle",[0 0 1 1],"Color",'w');
exportgraphics(fig,"ERA3SM-intercept.eps")
delete(a)


% Amplifier BFU-760F
L_att = [25 22 19 16 13 11 10 9 7 4];
P_a = [-4.5 -1.75 1.3 3.8 6 7 7.2 7.3 7.3 7.1];
P_IM3 = [-57 -48 -34.7 -23.6 -16.6 -16.6 -16.2 -16.1 -11.7 -5.7];
% O_IM3 = P_a - P_IM3; disp(O_IM3)

% approximate linear behvaiour at the beginning
coeff_a = polyfit(-L_att(1:5),P_a(1:5),1);
line_a = polyval(coeff_a,-L_att);
coeff_IM3 = polyfit(-L_att(1:5),P_IM3(1:5),1);
line_IM3 = polyval(coeff_IM3,-L_att);

% find IP3
X_alt = linspace(-5,5);
line_a_alt = coeff_a(1)*X_alt + coeff_a(2);
line_IM3_alt = coeff_IM3(1)*X_alt + coeff_IM3(2);
diff = line_a_alt - line_IM3_alt;
[~, i] = min(abs(diff));
IP3_in = round(X_alt(i),1);
IP3_out = round(line_a_alt(i),1);
disp(['LNA with BFU-760F: IP3 = ' num2str(IP3_out)])

% plot
fig = figure(4);
plot(-L_att,P_a,-L_att,P_IM3)
axis('tight')
hold on
plot(-L_att,line_a,'--',-L_att,line_IM3,'--','Color',[0.5 0.5 0.5])
% xline(IP3_in,'--','Color',[0.5 0.5 0.5])
% yline(IP3_out,'--','Color',[0.5 0.5 0.5])
grid on
grid minor
legend(["P_a" "P_{2a \pm b}"],"Location","Southeast")
% legend("boxoff")
xlabel("P_{in} [dBm]");
ylabel("P_{out} [dBm]");
% title("BFU-760F")
hold off

% export
% draw a white rectangle around the graph to avoid trimming during export
a = annotation("rectangle",[0 0 1 1],"Color",'w');
exportgraphics(fig,"BFU760F-intercept.eps")
delete(a)

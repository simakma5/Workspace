close all; clear; clc

%% Task 1: transform to LLH and plot estimated velocity in kmph
load('position_ecef.mat');
ECEFpos = xyz_record(:,2:1:4)';
llhpos = ecef2llh(ECEFpos,'wgs84')';

speed = sqrt( ...
    (xyz_record(2:end,2) - xyz_record(1:end-1,2)).^2 + ...
    (xyz_record(2:end,3) - xyz_record(1:end-1,3)).^2 + ...
    (xyz_record(2:end,4) - xyz_record(1:end-1,4)).^2)./abs( ...
    (xyz_record(2:end,1) - xyz_record(1:end-1,1)));

% speed = speed*3.6;
speed_avg = mean(speed);

windowSize = 12;
[b,a] = butter(2,0.03);
speed_filter = filter(b,a,speed);

figure(1)
plot(xyz_record(2:end,1),speed,"Color",[0.7 0.7 0.8])
hold on
plot(xyz_record(2:end,1),speed_filter,"LineWidth",0.7)
axis tight
xlabel("Time [h]")
ylabel("Velocity [km/h]")
title(['Average speed: ' num2str(speed_avg) ' km/h'])
hold off

% figure(11)
% plot(llhpos(:,3))
% hold on
% yyaxis right
% plot([0; speed_filter])
% hold off

%% Task 2: determine where the record originates from (country, city)
figure(2)
geoscatter(llhpos(:,1),llhpos(:,2),10,[0; speed_filter],'filled')
colormap("jet")

%% Task 3: display record in local coordinates
[xEast,yNorth,zUp] = geodetic2enu(llhpos(:,1),llhpos(:,2),llhpos(:,3),llhpos(1,1),llhpos(1,2),llhpos(1,3),wgs84Ellipsoid);
figure(3)
plot3(xEast, yNorth, zUp);
xlabel('x')
ylabel('y')
zlabel('z')

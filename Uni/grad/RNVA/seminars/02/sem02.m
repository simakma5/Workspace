close all; clear; clc
% 1: transform to LLH -> priblizny graf rychlosti v km/h
% 2: kde zaznam vznikl (zeme, mesto)
% 3: zobrazit v lok. sour. (xr,yr,zr); 1. bod zaznamu (0,0,0)_ENU
% bonus: zobrazeni vyskoveho profilu trati
%% Kubajz
load('position_ecef.mat');
ECEFpos = xyz_record(:,2:1:4);

llhpos=ecef2llh(ECEFpos','wgs84')';

windowSize = 12; 
[b, a] = butter(2, 0.03);

%writematrix(llhpos,'M.csv')
% figure(1)
% geo = geoplot(llhpos(:,1), llhpos(:,2));

wgs84 = wgs84Ellipsoid;

speed = sqrt((xyz_record(1:end-1,2) - xyz_record(2:end,2)).^2+(xyz_record(1:end-1,3) - xyz_record(2:end,3)).^2+(xyz_record(1:end-1,4) - xyz_record(2:end,4)).^2)./abs((xyz_record(1:end-1,1) - xyz_record(2:end,1)));

speed = speed*3.6;
avgspeed = mean(speed)

y = filter(b,a,speed);

figure(1)
plot(speed, "Color",[0.7 0.7 0.8])
hold on
plot(y,"LineWidth",0.7)
hold off

[xEast,yNorth,zUp] = geodetic2enu(llhpos(:,1), llhpos(:,2), llhpos(:,3), llhpos(1,1), llhpos(1,2), llhpos(1,3),wgs84);


figure(2)
geoscatter(llhpos(:,1), llhpos(:,2), 10, [0; y],'filled')
colormap("jet")

figure(3)
plot(llhpos(:,3))
hold on
yyaxis right
plot([0; y])
hold off

figure(4)
plot3(xEast, yNorth, zUp);
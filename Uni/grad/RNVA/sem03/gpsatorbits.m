%EASY11    Creates a stereographic plot of GPS satellite
%          orbits from an almanac. The plot is as seen
%          from the position (phi,lambda). All orbits
%          with elevation angle lower than a given
%          cut-off angle (mask) are omitted.
%          Almanac files most easily can be downloaded
%          from your own (high-end) receiver or
%          http://www.ngs.noaa.gov/CORS/Data.thml
%
%          An additional plot is created showing number
%          of visible satellites and when they are visible

% (PP)    Another plot shows trajectory of satellites in 3D ECEF coordinate
%         system vith camera position related to phi lambda reference

% based on original script easy11.m from
%Kai Borre 26-08-2008
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2008/08/26  $

%adjustments made by Pavel Puricer
%$Rev: 1.0.1 $ $Date: 2017/03/08  $

clear all
close all

set(0,'DefaultTextFontName','Times');
set(0,'DefaultAxesFontName','Times');
set(0,'DefaultTextFontSize',16);

% the following five assignments are necessary
% rinexe('6598067mr.17n','fel.nav');
% rinexe('log0307l_red.19N','fel.nav');
% almanac = 'fel.nav';
% elevation mask to exclude low elevation satellites
% because their information is likely damaged due to big troposcatter
mask = 5;
% Prague
phi = [50 6 0];
lambda = [14 23 0];
% Equator
% phi = [0 0 0];
% lambda = [0 0 0];
% North pole
% phi = [90 0 0];
% lambda = [0 0 0];
duration_sec=43200;
timestep_sec=600;
start_time=0;
gpsweek=1019;
leapsec=37;
origin_of_week=datenum(2023,04,16,0,0,0);
frac_of_week=start_time/86400;
timetag_of_start=origin_of_week+frac_of_week;
start_string=datestr(timetag_of_start);

%reading ephemerides
% fide = fopen(almanac,'r');
% Eph = fread(fide,inf,'double');
% m = length(Eph);
% eph = reshape(Eph,21,m/21);
% load('ephemerides_20170308.mat','eph')
eph = rinexe3("C:\Workspace\Uni\grad\RNVA\sem03\calculations\RINEX\GOP600CZE_R_20231100600_01H_GN.rnx");

% transformation of given location (phi,lambda,h) to (X,Y,Z)
Phi = dms2rad(phi(1),phi(2),phi(3));
Phi = Phi*180/pi;
Lambda = dms2rad(lambda(1),lambda(2),lambda(3));
Lambda = Lambda*180/pi;
[M(1,1),M(2,1),M(3,1)] = frgeod(6378137,298.257222101,Phi,Lambda,0);

% Computation of (azimuth, elevation) for each satellite
% for each 15 min.s. We use only one ephemeris for each PRN.
% Anyway, the plot is only for visual use
[prns, ind] = unique(eph(1,:));
az = ones(32,96)*inf;
el = ones(32,96)*inf;

for sat = [ind]'
    if isempty(start_time)
        start_time = eph(18,sat);
    end
    j = 0;
    i = eph(1,sat);
    for time = start_time:timestep_sec:start_time+duration_sec
        S = satpos(time,eph(:,sat));
        j = j+1;
        [azimuth,elevation,distance] = topocent(M,S-M);
        az(i,j) = azimuth;
        el(i,j) = elevation;
        possat_ecef(i,j,:)=S;
    end
end
%%
%%%%%%%%% stereographic plot  %%%%%%%%%%%%%
figure(1);
% polarhg draws coordinate lines of a polar plot. We add
% circles with radii 30 and 60 degrees
polarhg([30 60])
XX = zeros(32,40)*inf; % a matrix of NaNs to store plot data
YY = XX;

hold on
for k = 1:32
    if az(k,1) == 0, continue, end
    AZ = az(k,:);
    EL = el(k,:);
    % remove data below the cut-off angle
    AZ(find(EL <= mask)) = nan; 
    EL(find(EL <= mask)) = nan;
    % convertion to polar coordinates
    xx = (90-EL).*cos(AZ*pi/180);
    yy = (90-EL).*sin(AZ*pi/180);
    XX(k,1:length(xx)) = xx;
    YY(k,1:length(yy)) = yy;
end % end-k
% the first coord. moves text vertically (increasing values up),
% the second coord. moves text horizontally (increasing values right)
text(135,-95,{['Skyplot for the position (\phi, \lambda) = (' ...
    num2str(round(Phi)) '\circ, '  num2str(round(Lambda)) '\circ)']})
text(115,-45,{['Elevation mask  ' num2str(mask) '\circ' ]}) %120
text(-120,-120,['All PRNs except  ' num2str(setdiff(1:32,prns)) ])
plot(XX',YY','linewidth',2)

hold off

print -depsc2 easy111
%% Plot 3D view of satellite paths in ECEF
sphere_radius=6371000.7;
[spx,spy,spz]=sphere;
Espx=sphere_radius*spx;
Espy=sphere_radius*spy;
Espz=sphere_radius*spz;
figure(2)
hold on
mesh(Espx,Espy,Espz,ones(size(Espx)))
% uprava kvuli mensimu poctu druzic
dimension_ecef=size(possat_ecef);
%
for k=1:dimension_ecef(1,1)
    plot3(possat_ecef(k,:,1),possat_ecef(k,:,2),possat_ecef(k,:,3))
end
plot3([0 0],[-3e7 3e7],[0 0],'k')
plot3([-3e7 3e7],[0 0],[0 0],'k')
plot3([0 0],[0 0],[-3e7 3e7],'k')
hold off
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')

title('3D Satellites orbits in ECEF')
%plot3(possat_ecef(1,:,1),possat_ecef(1,:,2),possat_ecef(1,:,3))
campos(M)
print -depsc2 easy1110




%%
% preparation for visibility plot  %%%%%%%%%%%%%%%%%

% we choose a resolution of 5 min.s,
% ie. 24 hours times 12 = 288 which becomes the range of j
satsum = zeros(1,288);
visible = zeros(2*(size(prns,2)+1),288);

for sat = [ind]'
    Az = [];
    El = [];
    i = eph(1,sat);
    for j = 1:288
        time = 300*(j-1);
        S = satpos(start_time+time,eph(:,sat));
        [az,el,d] = topocent(M,S-M);
        if el > mask
            Az = [Az az];
            El = [El el];
            satsum(j) = satsum(j)+1;
            visible(2*i,j) = 1;
        end
    end
end

figure(3);
ax1=subplot(2,1,1)
set(gca,'Fontsize',16);
area(satsum)
colormap(ax1,gray)
set(gca,'XTick',1:71:288)
set(gca,'XTickLabel',{'0','6','12','18','24'})
xlabel('GPS Time [hours]')
ylabel('# of Visible Satellites')
title(['OriginTime: ',num2str(start_string),' Elevation Mask ' num2str(mask) '\circ'])

print -depsc2 easy112
%%
%figure(3);
ax2=subplot(2,1,2)
set(gca,'Fontsize',16);
imagesc(flipud(visible));
colormap(ax2,copper)
set(gca,'XTick',1:71:288)
set(gca,'XTickLabel',{'0','6','12','18','24'})
set(gca,'YTick',-3:16:(2*(size(prns,2)+1)))
set(gca,'YTickLabel',{'2','8','16','24','32'});
xlabel('GPS Time [hours]')
ylabel('PRNs')
title('Solid Lines Indicate Visible Satellites')
colormap summer

print -depsc2 easy113
%%%%%%%%%%%%%%%%%%%%% end easy11.m %%%%%%%%%%%%%%

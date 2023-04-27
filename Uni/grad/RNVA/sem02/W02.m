% ECEFpos = [0.01;0;6356752.31]; % North pole
% ECEFpos = [6378137;0;0]; % Greenwich+equator
ECEFpos = [3970599.791;1018943.437;4870317.63]; % FEL
disp(['Original input -  X: ',num2str(ECEFpos(1,1)),' Y: ',num2str(ECEFpos(2,1)),' Z: ',num2str(ECEFpos(3,1))])

llhpos=ecef2llh(ECEFpos,'wgs84');
disp(['Convereted - Lat: ',num2str(llhpos(1,1)),' Lon: ',num2str(llhpos(2,1)),' Alt: ',num2str(llhpos(3,1))])

ECEFpos2=llh2ecef(llhpos,'wgs84');
% ECEFpos2=llh2ecef(llhpos,'bessel1841');
% ECEFpos2=llh2ecef(llhpos,'krasov');
disp(['Backward conversion - X2: ',num2str(ECEFpos2(1,1)),' Y2: ',num2str(ECEFpos2(2,1)),' Z2: ',num2str(ECEFpos2(3,1))])


diff_x=ECEFpos(1,1)-ECEFpos2(1,1);
diff_y=ECEFpos(2,1)-ECEFpos2(2,1);
diff_z=ECEFpos(3,1)-ECEFpos2(3,1);

disp(['Error of fwd+back conversion [m] dx: ',num2str(diff_x),' dy: ',num2str(diff_y),' dz: ',num2str(diff_z)])

diff3d=sqrt(diff_x^2+diff_y^2+diff_z^2);
disp(['3D rms error of conversion forth&back [m]: ',num2str(diff3d)])

%% Porovnani rozdilu urceni polohy pri pouziti rozdilnych elisoidu (WGS84 vs. Krasovskij)

% llhpos10=[(0:0.1:90);zeros(1,length((0:0.1:90)));zeros(1,length((0:0.1:90)))];
lat=(0:0.1:90);
lon=45;
height=0;
llhpos10=[lat;lon*ones(1,length(lat));height*ones(1,length(lat))];
ECEFpos10wgs84=llh2ecef(llhpos10,'wgs84');
ECEFpos10nowgs=llh2ecef(llhpos10,'krasov');
% ECEFpos10nowgs=llh2ecef(llhpos10,'sphere');

figure(10)
subplot(3,1,1)
plot((0:0.1:90),ECEFpos10wgs84(1,:)-ECEFpos10nowgs(1,:))
title({'error of use of different ellipsoids:';['WGS 84 vs. Krasovskij Lon: ',num2str(lon),' Height: ',num2str(height)]})
ylabel('x [m]')
subplot(3,1,2)
plot((0:0.1:90),ECEFpos10wgs84(2,:)-ECEFpos10nowgs(2,:))
ylabel('y [m]')
subplot(3,1,3)
plot((0:0.1:90),ECEFpos10wgs84(3,:)-ECEFpos10nowgs(3,:))
ylabel('z [m]')
xlabel('Latitude [deg]')

figure(20)
diff3d_par=(ECEFpos10wgs84(1,:)-ECEFpos10nowgs(1,:)).^2+(ECEFpos10wgs84(2,:)-ECEFpos10nowgs(2,:)).^2+(ECEFpos10wgs84(3,:)-ECEFpos10nowgs(3,:)).^2;
plot((0:0.1:90),sqrt(diff3d_par))
title({'3d RMS error of use of different ellipsoids:';'WGS84 vs. Krasovskij'})
ylabel('3D RMS [m]')
xlabel('Latitude [deg]')

% figure(11)
% plot3(ECEFpos10wgs84(1,:),ECEFpos10wgs84(2,:),ECEFpos10wgs84(3,:))
% hold on
% plot3(ECEFpos10nowgs(1,:),ECEFpos10nowgs(2,:),ECEFpos10nowgs(3,:))
% hold off


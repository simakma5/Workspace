function satpos = satpos(t,eph)
% Computation of GPS and Galileo Coordinates according to
% https://gssc.esa.int/navipedia/index.php/GPS_and_Galileo_Satellite_Coordinates_Computation

% From rinexe.m:
% eph(1,:)  = svprn;
% eph(2,:)  = af2;
% eph(3,:)  = M0;
% eph(4,:)  = roota;
% eph(5,:)  = deltan;
% eph(6,:)  = ecc;
% eph(7,:)  = omega;
% eph(8,:)  = cuc;
% eph(9,:)  = cus;
% eph(10,:) = crc;
% eph(11,:) = crs;
% eph(12,:) = i0;
% eph(13,:) = idot;
% eph(14,:) = cic;
% eph(15,:) = cis;
% eph(16,:) = Omega0;
% eph(17,:) = Omegadot;
% eph(18,:) = toe;
% eph(19,:) = af0;
% eph(20,:) = af1;
t_oe = eph(18,:);               % ephemerides reference epoch
roota = eph(4,:);               % square root of semi-major axis
ecc = eph(6,:);                 % eccentricity
M0 = eph(3,:);                  % mean anomaly at reference epoch
omega = eph(7,:);               % argument of perigee
i0 = eph(12,:);                 % inclination at reference epoch
Omega0 = eph(16,:);
delta_n = eph(5,:);             % mean motion difference
i_dot = eph(13,:);              % rate of inclination angle
Omega_dot = eph(17,:);          % rate of node's right ascension
cuc = eph(8,:);                 % latitude argument correction
cus = eph(9,:);                 % latitude argument correction
crc = eph(10,:);                % orbital radius correction
crs = eph(11,:);                % orbital radius correction
cic = eph(14,:);                % inclination correction
cis = eph(15,:);                % inclination correction

mu = 3.986005e14;               % gravitational parameter of Earth
omegaE = 7.2921151467e-5;       % Earth's rotation

tk = check_t(t-t_oe);
n0 = sqrt(mu/roota^6);
Mk = M0+(n0+delta_n)*tk;
Mk = rem(Mk+2*pi,2*pi);
Ek = Mk;
for i = 1:10
    E_old = Ek;
    Ek = Mk+ecc*sin(Ek);
    if abs(rem(Ek-E_old,2*pi)) < 1e-12
        break;
    end
end
Ek = rem(Ek+2*pi,2*pi);
vk = atan2(sqrt(1-ecc^2)*sin(Ek),cos(Ek)-ecc);
phi = omega+vk;
phi = rem(phi,2*pi);
uk = phi+cuc*cos(2*phi)+cus*sin(2*phi);
rk = roota^2*(1-ecc*cos(Ek))+crc*cos(2*phi)+crs*sin(2*phi);
ik = i0+i_dot*tk+cic*cos(2*phi)+cis*sin(2*phi);
lambdak = Omega0+(Omega_dot-omegaE)*tk-omegaE*t_oe;
lambdak = rem(lambdak+2*pi,2*pi);

satpos(1,1) = rk*(cos(lambdak)*cos(uk)-cos(ik)*sin(lambdak)*sin(uk));
satpos(2,1) = rk*(sin(lambdak)*cos(uk)+cos(ik)*cos(lambdak)*sin(uk));
satpos(3,1) = rk*sin(ik)*sin(uk);

% x1 = cos(uk)*rk;
% y1 = sin(uk)*rk;
% satpos(1,1) = x1*cos(lambdak)-y1*cos(ik)*cos(lambdak);
% satpos(2,1) = x1*sin(lambdak)+y1*cos(ik)*cos(lambdak);
% satpos(3,1) = y1*sin(ik);

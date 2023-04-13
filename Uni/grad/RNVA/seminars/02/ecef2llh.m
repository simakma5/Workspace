function llh = ecef2lla(ecefpos,ellips)
%
% ecef2lla
% Converts ecef coordinates to lla coordinates.
% [latitude, longitude, altitude] = ecef2lla(ECEF position)
% ECEF Position Array = [x; y; z]
% x, y, z are row vectors with same length
% 
% 
% 
% A negative longitude is the Western hemispere.
% A negative latitude is in the Southern hemisphere.
% 
% Ellipsoids are loaded from exterenal file Ellipsoids2.mat
% WGS84 parameters
% a = 6378137.0; %semimajor axis (equatorial) radius
% b = 6356752.3142; % semiminor axis (polar) radius

%% Load ellipsoids
load Ellipsoids2;
if ~exist(ellips,'var')
    error(['Ellipsoid ',ellips,' is not defined in Ellipsoids.mat - check your definitions!.'])
end
eval(['ell=',ellips,';']);
a=ell.a;
b=ell.b;
% f=ell.f;
%%

f = (a-b)/a; %flattenning
e = sqrt(f*(2-f)); % eccentricity of ellipsoid
len = length(ecefpos(1,:)); %get length of data
for i = 1:len
    lon(i) = atan2(ecefpos(2,i), ecefpos(1,i));
    lon(i) = lon(i)*180/pi; %convert to degrees

    h = 0; % initialization
    N = a;
    flag = 0;
    j = 0;
    p = sqrt(ecefpos(1,i)^2 + ecefpos(2,i)^2);
    
    sinlat = ecefpos(3,i)/(N*(1-e^2)+h);  %First iteration
    lat(i) = atan((ecefpos(3,i)+e^2*N*sinlat)/p);
    N = a/(sqrt(1 - (e^2)*(sinlat^2)));
    prevalt = (p/cos(lat(i)))-N;
    prevlat = lat(i)*180/pi;
    
    while (flag < 2) %do at least 100 iterations
        flag = 0;
        sinlat = ecefpos(3,i)/(N*(1-e^2)+h);
        lat(i) = atan((ecefpos(3,i)+e^2*N*sinlat)/p);
        N = a/(sqrt(1 - (e^2)*(sinlat^2)));
        alt(i) = (p/cos(lat(i)))-N;
        lat(i) = lat(i)*180/pi;
        if abs(prevalt-alt(i)) < .00000001
            flag = 1;
        end
        if abs(prevlat-lat(i)) < .00000001
            flag = flag + 1;
        end
        j = j+1;
        if j == 100
            flag = 2;
        end
        prevalt = alt(i);
        prevlat = lat(i);
    end
end

llh=[lat;lon;alt];

return


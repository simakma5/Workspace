function ecf=llh2ecef(llh,ellips)
%
% Description:
% 		This function returns the x, y, and z earth centered fixed (ECF)
%		coordinates of the point specified by llh [lat lon hgt].  
%		Note that longitude is positive east of the Greenwich meridian.
%
% Usage:
%   ecf=llh2ecef(llh,ellips)
%
% Inputs:
%   llh   - 3xN Matrix of Vectors
%           where
%              llh(1,:) is latitude positive (degrees)
%              llh(2,:) is longitude positive East (degrees)
%              llh(3,:) is height (meters)
%   elips  - string with ellipsoid name, e.g. 'wgs84','krasov','bessel1841'
%
% Outputs
%   ecf   - 3xN Matrix of Vectors in ECF coordinates (meters)
%           where
%              x(1,:) is Greenwich meridan; (0 lon, 0 lat)
%              y(2,:) is Easterly
%              z(3,:) is North Pole
%

%% Load ellipsoids
load Ellipsoids2;
if ~exist(ellips,'var')
    error(['Ellipsoid ',ellips,' is not defined in Ellipsoids2.mat - check your definitions!.'])
end
eval(['ell=',ellips,';']);
a=ell.a;
b=ell.b;
f=ell.f;
%%


if isvector(llh)
    llh = llh(:);
end

lat=llh(1,:)*(pi/180);
lon=llh(2,:)*(pi/180);
hgt=llh(3,:);

%  Set up WGS-84 constants.
%[a,f] = EarthModel_sphere;

%  Convert lat,lon,hgt geographic coords to XYZ Earth Centered Earth Fixed coords.
%		N = a/sqrt( 1 - f*(2-f)*sin(lat)*sin(lat) )
%		X = (N + h)*cos(lat)*cos(lon)
%		Y = (N + h)*cos(lat)*sin(lon)
%		Z = ((1-f)^2 * N + h)*sin(lat)

%  Store some commonly used values.

slat = sin(lat);
N = a./sqrt(1 - f*(2-f) * slat.^2);
Nplushgtclat = (N + hgt) .* cos(lat);

x = Nplushgtclat .* cos(lon);
y = Nplushgtclat .* sin(lon);
z = ((1-f)^2 * N + hgt) .* slat;

ecf = [x; y; z];

return

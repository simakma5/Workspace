function eph = rinexe3(ephemerisfile)
%RINEXE Reads a RINEX Navigation Message file version 3 and
%	    reformats the data into a matrix eph with 21
%	    rows and a column for each satellite.


%Typical call: ephem = rinexe('pta.96n')

%Kai Borre 04-18-96
%Copyright (c) by Kai Borre
%$Revision: 1.1 $  $Date: 2008/06/02  $

%Modification Pavel Puricer 2023-04-20
% input: RINEX nav file
% output: variable eph with ephemerides

% Units are either seconds, meters, or radians
fide = fopen(ephemerisfile);
head_lines = 0;
while 1  % We skip header
   head_lines = head_lines+1;
   line = fgetl(fide);
   answer = findstr(line,'END OF HEADER');
   if ~isempty(answer), break;	end;
end;
disp(['head_lines ', num2str(head_lines)])
noeph = -1;
while 1
   noeph = noeph+1;
   line = fgetl(fide);
   if line == -1, break;  end
end;
noeph = noeph/8;
disp(['No. of ephemerides: ',num2str(noeph)])
frewind(fide);
for i = 1:head_lines, line = fgetl(fide); end;

% Set aside memory for the input
svprn	 = zeros(1,noeph);
weekno	 = zeros(1,noeph);
t0c	 = zeros(1,noeph);
tgd	 = zeros(1,noeph);
aodc	 = zeros(1,noeph);
toe	 = zeros(1,noeph);
af2	 = zeros(1,noeph);
af1	 = zeros(1,noeph);
af0	 = zeros(1,noeph);
aode	 = zeros(1,noeph);
deltan	 = zeros(1,noeph);
M0	 = zeros(1,noeph);
ecc	 = zeros(1,noeph);
roota	 = zeros(1,noeph);
toe	 = zeros(1,noeph);
cic	 = zeros(1,noeph);
crc	 = zeros(1,noeph);
cis	 = zeros(1,noeph);
crs	 = zeros(1,noeph);
cuc	 = zeros(1,noeph);
cus	 = zeros(1,noeph);
Omega0	 = zeros(1,noeph);
omega	 = zeros(1,noeph);
i0	 = zeros(1,noeph);
Omegadot = zeros(1,noeph);
idot	 = zeros(1,noeph);
accuracy = zeros(1,noeph);
health	 = zeros(1,noeph);
fit	 = zeros(1,noeph);

for i = 1:noeph
   line = fgetl(fide);	  %
   svprn(i) = str2num(line(2:3));
   year(i) = str2num(line(5:8));
   month(i) = str2num(line(10:11));
   day(i) = str2num(line(13:14));
   hour(i) = str2num(line(16:17));
   minute(i) = str2num(line(19:20));
   second(i) = str2num(line(22:23));
   af0(i) = str2num(line(24:42));
   af1(i) = str2num(line(43:61));
   af2(i) = str2num(line(62:80));
   line = fgetl(fide);	  %
   IODE(i) = str2num(line(5:23));
   crs(i) = str2num(line(24:42));
   deltan(i) = str2num(line(43:61));
   M0(i) = str2num(line(62:80));
   line = fgetl(fide);	  %
   cuc(i) = str2num(line(5:23));
   ecc(i) = str2num(line(24:42));
   cus(i) = str2num(line(43:61));
   roota(i) = str2num(line(62:80));
   line=fgetl(fide);
   toe(i) = str2num(line(5:23));
   cic(i) = str2num(line(24:42));
   Omega0(i) = str2num(line(43:61));
   cis(i) = str2num(line(62:80));
   line = fgetl(fide);	    %
   i0(i) =  str2num(line(5:23));
   crc(i) = str2num(line(24:42));
   omega(i) = str2num(line(43:61));
   Omegadot(i) = str2num(line(62:80));
   line = fgetl(fide);	    %
   idot(i) = str2num(line(5:23));
   codes = str2num(line(24:42));
   weekno = str2num(line(43:61));
   L2flag = str2num(line(62:80));
   line = fgetl(fide);	    %
   svaccur = str2num(line(5:23));
   svhealth = str2num(line(24:42));
   tgd(i) = str2num(line(43:61));
   iodc = line(62:80);
   line = fgetl(fide);	    %
   tom(i) = str2num(line(5:23));
%   spare = line(24:42); % recent RINEX files have empty spaces
%   spare = line(43:61);
%   spare = line(62:80);
end
status = fclose(fide);

%  Description of variable eph.
eph(1,:)  = svprn; % satellite PRN number
eph(2,:)  = af2; % SV clock drift rate
eph(3,:)  = M0; % Mean anomaly at reference epoch
eph(4,:)  = roota; % sqrt of semi-major axis
eph(5,:)  = deltan; % mean motion difference
eph(6,:)  = ecc; % eccentricity
eph(7,:)  = omega; % Rate of node's right ascension
eph(8,:)  = cuc; % Latitude argument correction cuc
eph(9,:)  = cus; % Latitude argument correction cus
eph(10,:) = crc; % Orbital radius correction crc
eph(11,:) = crs; % Orbital radius correction crs
eph(12,:) = i0; % Inclination at reference epoch
eph(13,:) = idot; % Rate of inclination angle
eph(14,:) = cic; % Inclination correction cic
eph(15,:) = cis; % Inclination correction cis
eph(16,:) = Omega0; % Longitude of ascending node at the beginning of the week
eph(17,:) = Omegadot; % Rate of node's right ascension
eph(18,:) = toe; % Eph. reference epoch in seconds within the week
eph(19,:) = af0; % SV clock offset
eph(20,:) = af1; % SV clock drift
eph(21,:) = year*1e10+month*1e8+day*1e6+hour*1e4+minute*1e2+second; % TOC - Time of clock, epoch YYMMDDhhmmss.s

% fidu = fopen(outputfile,'w');
% count = fwrite(fidu,[eph],'double');
fclose all;
disp([num2str(length(eph(1,:))), ' Satellites loaded: ',num2str(eph(1,:))])
%%%%%%%%% end rinexe.m %%%%%%%%%
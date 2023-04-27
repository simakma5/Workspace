%elevace a azimut pouzitych druzic
%ea = [60,067;23,103;25,057;18,168;30,104;12,324];
ea = [5,0;5,120;5,240;90,0];
%ea = [0,0;0,90;0,180;1,270];
%ea = [0,0;0,90;0,180;45,270;70,100];
%ea = [90,0;-30,0;-30,120;-30,240];
% $GPGSA,A,3,21,06,24,22,07,25,,,,,,,2.5,1.2,2.1*36
% $GPGSV,3,1,12,16,67,237,,18,40,131,,03,35,288,,21,60,067,41*70
% $GPGSV,3,2,12,06,23,103,35,19,07,282,22,24,25,057,42,22,18,168,26*78
% $GPGSV,3,3,12,07,30,104,38,27,05,342,30,25,12,324,34,26,04,064,39*78
% ea = [60,067;23,103;25,057;18,168;30,104;12,324];
% vysledek dops =    2.7280    2.5479    1.2509    2.2197

% $GPGSA,A,1,06,07,10,24,,,,,,,,,26.4,18.1,19.2*04
% $GPGSV,2,1,07,13,05,346,26,12,01,128,,06,69,071,36,16,34,302,*7D
% $GPGSV,2,2,07,07,77,064,38,10,21,058,37,24,64,097,37*45
% ea = [69,071;77,064;21,058;64,097];
% vysledk dops =   35.6532   26.8141   18.3453   19.5562   23.4979

% $GPGSA,A,3,21,06,07,10,24,,,,,,,,7.2,4.8,5.3*38
% $GPGSV,2,1,05,21,75,114,33,06,42,091,38,07,50,091,37,10,12,032,34*73
% $GPGSV,2,2,05,24,45,059,42*41
% ea = [75,114;42,091;50,091;12,032;45,059];
% vysledek dops =   8.9601    7.1284    4.8070    5.2638    5.4285

ea = ea/180*pi; %prevod na radiany
los = [sin(ea(:,2)).*cos(ea(:,1)),cos(ea(:,2)).*cos(ea(:,1)),sin(ea(:,1))]; %jednotkove vektory line-of-sight
c = [los ones(size(los,1),1)]

no_sats = length(ea);
X = zeros(4,1);
Y = zeros(4,1);
Z = zeros(4,1);
figure(1)
quiver3(X,Y,Z,los(:,1),los(:,2),los(:,3),'*')
figure(2)
polar(ea(:,2),cos(ea(:,1)),'*')

dop = inv(c'* c);

dops(1) = sqrt(dop(1,1) + dop(2,2) + dop(3,3) + dop(4,4));    %  GDOP
dops(2) = sqrt(dop(1,1) + dop(2,2) + dop(3,3));               %  PDOP
dops(3) = sqrt(dop(1,1) + dop(2,2));                          %  HDOP
dops(4) = sqrt(dop(3,3));                                     %  VDOP 
dops(5) = sqrt(dop(4,4));                                     %  TDOP

disp(['GDOP  ','PDOP  ','HDOP  ','VDOP  ','TDOP'])
disp([num2str(dops(1),3),'  ',num2str(dops(2),3),'  ',num2str(dops(3),3),'  ',num2str(dops(4),3),'  ',num2str(dops(5),3)])
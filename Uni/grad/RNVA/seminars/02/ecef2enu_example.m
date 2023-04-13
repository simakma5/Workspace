% llh_origin = [14.39278, 50.10278, 245]; % FEL LLH
% xyz_origin = [3970599.791,1018946.437,4870317.63]; %FEL xyz (WGS84)
llh_origin = [0, 0, 0]; % Greenwich zero meridian, equator
xyz_origin = [6378137,0,0];

lambda_o = llh_origin(1,1);
phi_o = llh_origin(1,2);
trans_ecef2enu = [-1*sin(pi/180*lambda_o), cos(pi/180*lambda_o),0;-1*cos(pi/180*lambda_o)*sin(pi/180*phi_o), -1*sin(pi/180*lambda_o)*sin(pi/180*phi_o), cos(pi/180*phi_o);cos(pi/180*lambda_o)*cos(pi/180*phi_o),sin(pi/180*lambda_o)*cos(pi/180*phi_o), sin(pi/180*phi_o)];
% inversion using simple transposition
xyz_src = [3970599.791,1018946.437,4870317.63];
% xyz_delta = xyz_src - xyz_origin;
xyz_delta = [100, 100, 0];
ENU = (trans_ecef2enu*xyz_delta')';

disp(['local origin in LLH: Lat: ',num2str(llh_origin(1,1)),' Lon: ',num2str(llh_origin(1,2)),' HAE: ',num2str(llh_origin(1,3))])
disp(['local origin in ECEF: x: ',num2str(xyz_origin(1,1)),' y: ',num2str(xyz_origin(1,2)),' z: ',num2str(xyz_origin(1,3))])
disp(['ECEF relative point: x: ',num2str(xyz_delta(1,1)),' y: ',num2str(xyz_delta(1,2)),' z: ',num2str(xyz_delta(1,3))])
disp(['ENU position: East: ',num2str(ENU(1,1)),' North: ',num2str(ENU(1,2)),' Up: ',num2str(ENU(1,3))])

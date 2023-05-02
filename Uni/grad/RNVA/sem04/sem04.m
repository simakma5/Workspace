close all; clear; clc
load("C:/Workspace/Uni/grad/RNVA/sem04/seminarw10_data.mat")
% cas_sec       time within GPS week [s]
% p_ranges      pseudodistances [m]
% eph           satellites ephemerides

% cas_sec,eph -> [x_ik,y_ik,z_ik] (i-th satellite, time t_k)
% p_ranges,[x_ik,y_ik,z_ik] -> x,y,z (user)

%% 1: using satpos(), determine x,y,z of satellites for given time t_k
% watch out: i != eph(1,i)
% this means that the index doesn't correspond to the satellite number
xyz_sat_1 = satpos(cas_sec(1,1),eph(:,1));
xyz_sat_5 = satpos(cas_sec(1,1),eph(:,2));
xyz_sat_10 = satpos(cas_sec(1,1),eph(:,3));
xyz_sat_12 = satpos(cas_sec(1,1),eph(:,4));

%% 2: choose suitable satellites from previous step (at least 4) and corresponding pseudodistances
D1 = p_ranges(1,1);     % satellite PRN1
D2 = p_ranges(1,5);     % satellite PRN5
D3 = p_ranges(1,10);    % satellite PRN10
D4 = p_ranges(1,12);    % satellite PRN12

%% 3: apply the iterational process for determination of position x,y,z (ECEF)
rv1 = xyz_sat_1;
rv2 = xyz_sat_5;
rv3 = xyz_sat_10;
rv4 = xyz_sat_12;

r0 = [0;0;0];
step_modifier = 0.1;
% step_modifier = 0.5;
% step_modifier = 1;
iter_steps = round(20/step_modifier);

for n=1:iter_steps
    d_od1 = sqrt((r0-rv1)'*(r0-rv1));
    d_od2 = sqrt((r0-rv2)'*(r0-rv2));
    d_od3 = sqrt((r0-rv3)'*(r0-rv3));
    d_od4 = sqrt((r0-rv4)'*(r0-rv4));
    A = [(r0-rv1)'/d_od1,1;(r0-rv2)'/d_od2,1;(r0-rv3)'/d_od3,1;(r0-rv4)'/d_od4,1];
    B = [D1-d_od1;D2-d_od2;D3-d_od3;D4-d_od4];
    rx = inv(A)*inv(A')*A'*B;
    r0 = r0+step_modifier*rx(1:3);
    disp(['Iteration: ', num2str(n), ' out of ', num2str(iter_steps)])
    disp(['Estimated r0: ', num2str(r0(1)),',',num2str(r0(2)),',',num2str(r0(3))])
    disp(['Step rx: ', num2str(rx(1)), ', ', num2str(rx(2)), ', ', num2str(rx(3))])
    disp('----')
    if rx(1) < 1e-1 && rx(2) < 1e-1 && rx(3) < 1e-1
        disp('Termination condition reached')
        break;
    end
end

%% 3a: EFEC to LLH -> geographic coordinates
llhpos = ecef2llh(r0,'wgs84');
figure(1)
geoscatter(llhpos(1),llhpos(2))
colormap("jet")

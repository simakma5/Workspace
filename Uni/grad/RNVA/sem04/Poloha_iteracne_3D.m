clear all
close all

%speed of light constant definition
c = 299792458;

%Polohy vysilacu v [m]
rv1 = [8000 8000 15000]';
% rv2 = [0 0 0]';
rv2 = [-2000 -2000 0]';
rv3 = [16000 0 0]';
rv4 = [8000 16000 0]';

rr = [4000 3500 800]'; %skutecna poloha
% rb = 0.06/1000*c; %skutecny posun zakladny

deltasys=6e-6; %posun casove zakladny [s]
%  delay_meas_err_max=1.1e-7; %chyba urceni zpozdeni [s]
delay_meas_err_max=0.11e-7; %chyba urceni zpozdeni [s]
% delay_meas_err_max=0.011e-7; %chyba urceni zpozdeni [s]
% delay_meas_err_max=0; %chyba urceni zpozdeni [s]

%Pseudovzdalenosti v m (puvodni verze, nzni pocitano na zaklade skutecnych
%vzdalenosti a chyby mereni pseudorange(delay))
% D1 = 0.11156/1000*c;
% D2 = 0.078222/1000*c;
% D3 = 0.102/1000*c;
% D4 = 0.104/1000*c;

% Matrix with summarized transmitter coordinates
Tx_matrix=[rv1,rv2,rv3,rv4];
% |x1 x2 x3 x4|
% |y1 y2 y3 y4|
% |z1 z2 z3 z4|

for (i=1:4)
    dist(i)=sqrt((Tx_matrix(:,i)-rr)'*(Tx_matrix(:,i)-rr));
    delay(i)=dist(i)/c;
    pseudodelay(i)=delay(i)+deltasys;
    delay_meas_err(i)=delay_meas_err_max*sign(rand-0.5);
    mpseudodelay(i)=pseudodelay(i)+delay_meas_err(i);
end

% dist1=sqrt((rv1-rr)'*(rv1-rr));
% delay1=dist1/c;
% pseudodelay1=delay1+deltasys;
% delay1_meas_err=delay_meas_err;
% mpseudodelay1=pseudodelay1+delay1_meas_err;
% 
% dist2=sqrt((rv2-rr)'*(rv2-rr));
% delay2=dist2/c;
% pseudodelay2=delay2+deltasys;
% delay2_meas_err=delay_meas_err;
% mpseudodelay2=pseudodelay2+delay2_meas_err;
% 
% dist3=sqrt((rv3-rr)'*(rv3-rr));
% delay3=dist3/c;
% pseudodelay3=delay3+deltasys;
% delay3_meas_err=delay_meas_err;
% mpseudodelay3=pseudodelay3+delay3_meas_err;
% 
% dist4=sqrt((rv4-rr)'*(rv4-rr));
% delay4=dist4/c;
% pseudodelay4=delay4+deltasys;
% delay4_meas_err=delay_meas_err;
% mpseudodelay4=pseudodelay4+delay4_meas_err;
% 
% D1 = mpseudodelay1*c;
% D2 = mpseudodelay2*c;
% D3 = mpseudodelay3*c;
% D4 = mpseudodelay4*c;

D1 = mpseudodelay(1)*c;
D2 = mpseudodelay(2)*c;
D3 = mpseudodelay(3)*c;
D4 = mpseudodelay(4)*c;


% D1 = 0.0925/1000*c;
% D2 = 0.1445/1000*c;
% D3 = 0.137/1000*c;
% D4 = 0.16525/1000*c;
% rr = [1000 3000 100]; %skutecna poloha
% rb = 0.08/1000*c; %skutecny posun zakladny


% r0 = [0;0;0];
% r0 = [10000;10000;5000];
r0 = [11111;12999;10111];
% r0 = [11111;17999;1111];
% r0 = [11111111;1111111;11111];
r_orig=r0;
r0_res=r0;

% step_modifier=0.1;
step_modifier=0.5;
% step_modifier=1;
% iter_steps=20;
iter_steps=round(20/step_modifier);

Err_pos_RMS=zeros(1,iter_steps);
step_rms=zeros(1,iter_steps);
for n=1:iter_steps
    d_od1 = sqrt((r0-rv1)'*(r0-rv1));
    d_od2 = sqrt((r0-rv2)'*(r0-rv2));
    d_od3 = sqrt((r0-rv3)'*(r0-rv3));
    d_od4 = sqrt((r0-rv4)'*(r0-rv4));
    A = [(r0-rv1)'/d_od1,1;(r0-rv2)'/d_od2,1;(r0-rv3)'/d_od3,1;(r0-rv4)'/d_od4,1];
    B = [D1-d_od1;D2-d_od2;D3-d_od3;D4-d_od4];
%     rx1 = A\B;
%     rx2 = inv(A)*B;
    rx = inv(A'*A)*A'*B;
    step_rms(n)=sqrt(sum(rx(1:3).^2));
%     rx-rx1
%     rx-rx2
    r0 = step_modifier*rx(1:3,:)+r0;
%     r0 = rx(1:3,:)/10+r0;
    r0_res(:,n)=r0;
    Err_pos_RMS(n)=sqrt((r0-rr)'*(r0-rr));
    disp(['Iteration: ',num2str(n),',estimated r0: ',num2str(r0(1)),',',num2str(r0(2)),',',num2str(r0(3)),' RMS of pos. error: ',num2str(Err_pos_RMS(n))])
    disp(['rx step: ',num2str(rx(1)),',',num2str(rx(2)),',',num2str(rx(3)),' step RMS: ',num2str(step_rms(n))])
end
%% Plotting the values of the variables during iterations
figure(1)
plot(log10(step_rms))
xlabel('Iteration count')
ylabel('log[m]')
title('RMS of iteration step')
%% Plot 3D image with stations, true position and calculated positions
figure(2)
hold on
plot3(rv1(1,1),rv1(2,1),rv1(3,1),'ob')
plot3(rv2(1,1),rv2(2,1),rv2(3,1),'og')
plot3(rv3(1,1),rv3(2,1),rv3(3,1),'om')
plot3(rv4(1,1),rv4(2,1),rv4(3,1),'ok')
plot3(rr(1,1),rr(2,1),rr(3,1),'or')
plot3(r_orig(1,1),r_orig(2,1),r_orig(3,1),'*r')
x_axis_limit=[min([rv1(1,1),rv2(1,1),rv3(1,1),rv4(1,1)]) max([rv1(1,1),rv2(1,1),rv3(1,1),rv4(1,1)])];
y_axis_limit=[min([rv1(2,1),rv2(2,1),rv3(2,1),rv4(2,1)]) max([rv1(2,1),rv2(2,1),rv3(2,1),rv4(2,1)])];
z_axis_limit=[min([rv1(3,1),rv2(3,1),rv3(1,1),rv4(3,1)]) max([rv1(3,1),rv2(3,1),rv3(3,1),rv4(3,1)])];
plot3(x_axis_limit,[0 0],[0 0],'k')
plot3([0 0],y_axis_limit,[0 0],'k')
plot3([0 0],[0 0],z_axis_limit,'k')

% plot3(r0_res(1,:),r0_res(2,:),r0_res(3,:),'x-')
plot3([r_orig(1,1),r0_res(1,:)],[r_orig(2,1),r0_res(2,:)],[r_orig(3,1),r0_res(3,:)],'x-')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
legend('TX1','TX2','TX3','TX4','RX','R0')
title('Position estimate')
hold off
campos([1e5,1e5,1e5])
figure(3)
subplot(2,2,3)
hold on
plot(r0_res(1,:),r0_res(2,:))
% plot(rv1(1,1),rv1(2,1),'ob')
% plot(rv2(1,1),rv2(2,1),'og')
% plot(rv3(1,1),rv3(2,1),'om')
% plot(rv4(1,1),rv4(2,1),'ok')
plot(rr(1,1),rr(2,1),'o')
xlabel('x [m]')
ylabel('y [m]')
hold off
title('Position in XY plane')
subplot(2,2,1)
hold on
plot(r0_res(1,:),r0_res(3,:))
% plot(rv1(1,1),rv1(3,1),'ob')
% plot(rv2(1,1),rv2(3,1),'og')
% plot(rv3(1,1),rv3(3,1),'om')
% plot(rv4(1,1),rv4(3,1),'ok')
plot(rr(1,1),rr(3,1),'o')
xlabel('x [m]')
ylabel('z [m]')
hold off
title('Position in XZ plane')
subplot(2,2,4)
hold on
plot(r0_res(3,:),r0_res(2,:))
% plot(rv1(3,1),rv1(2,1),'ob')
% plot(rv2(3,1),rv2(2,1),'og')
% plot(rv3(3,1),rv3(2,1),'om')
% plot(rv4(3,1),rv4(2,1),'ok')
plot(rr(3,1),rr(2,1),'o')
xlabel('z [m]')
ylabel('y [m]')
hold off
title('Position in ZY plane')

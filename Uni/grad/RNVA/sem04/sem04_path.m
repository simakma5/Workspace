close all; clear; clc
load("C:/Workspace/Uni/grad/RNVA/sem04/seminarw10_data.mat")

%%
pos = zeros(3,length(cas_sec));
r0 = [0;0;0];
for k=1:length(cas_sec)
    % satellite PRN1
    rv1 = satpos(cas_sec(1,k),eph(:,1));
    D1 = p_ranges(k,1);
    % satellite PRN5
    rv2 = satpos(cas_sec(1,k),eph(:,2));
    D2 = p_ranges(k,5);
    % satellite PRN10
    rv3 = satpos(cas_sec(1,k),eph(:,3));
    D3 = p_ranges(k,10);
    % satellite PRN12
    rv4 = satpos(cas_sec(1,k),eph(:,4));
    D4 = p_ranges(k,12);
    
    step_modifier=0.1;
%     step_modifier=0.5;
%     step_modifier=1;
%     iter_steps=20;
    iter_steps=round(20/step_modifier);
    
    for n=1:iter_steps
        d_od1 = sqrt((r0-rv1)'*(r0-rv1));
        d_od2 = sqrt((r0-rv2)'*(r0-rv2));
        d_od3 = sqrt((r0-rv3)'*(r0-rv3));
        d_od4 = sqrt((r0-rv4)'*(r0-rv4));
        A = [(r0-rv1)'/d_od1,1;(r0-rv2)'/d_od2,1;(r0-rv3)'/d_od3,1;(r0-rv4)'/d_od4,1];
        B = [D1-d_od1;D2-d_od2;D3-d_od3;D4-d_od4];
        rx = inv(A)*inv(A')*A'*B;
        r0 = r0+step_modifier*rx(1:3);
        disp(['Time index k: ', num2str(k), ' out of ', num2str(length(cas_sec))])
        disp(['Iteration: ', num2str(n), ' out of ', num2str(iter_steps)])
        disp(['Estimated r0: ', num2str(r0(1)),',',num2str(r0(2)),',',num2str(r0(3))])
        disp(['Step rx: ', num2str(rx(1)), ', ', num2str(rx(2)), ', ', num2str(rx(3))])
        disp('----')
        if rx(1) < 1e-1 && rx(2) < 1e-1 && rx(3) < 1e-1
            disp('Termination condition reached')
            break;
        end
    end
    pos(:,k) = r0;
end

%%
llhpos = ecef2llh(pos,'wgs84');
fig = figure(1);
geoscatter(llhpos(1,:),llhpos(2,:),10,'filled')
colormap("jet")
% draw a white rectangle around the graph to avoid trimming during export
a = annotation("rectangle",[0 0 1 1],"Color",'w');
exportgraphics(fig,'sem04_path.eps')
delete(a)

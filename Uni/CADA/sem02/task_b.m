%%
close all; clear; clc

%% initialization
% electrode potentials
V_0 = 100;
V_gnd = 0;
% permittivities
eps_r1 = 1;
eps_r2 = 4;
eps_av = (eps_r1+eps_r2)/2;
% grid
N_x = 20;
N_y = N_x;
N_if1 = round(3*N_x/4);
N_if2 = round(N_x/2);

% potential initialization
V = zeros(N_y,N_x);
V(1,:) = V_gnd;
V(N_x,:) = V_0;
V(:,1) = V_gnd;
V(:,N_y) = V_gnd; 
for i = 2:N_y-1
    for j = 2:N_x-1
        V(i,j) = 0;
    end
end
    
%% Solution of the trough problem using Successive over Relaxation (SOR) method
% update potencial values at free nodes
% SOR coefficient omega is choosen from the interval <1;2)
omega = 1.5;
% desired precision threshold
treshold = 0.001;
Rmax = V_0;
itnumber = 0;
R = zeros(N_y,N_x);
while Rmax > treshold
    for i = 2:N_y-1
        for j = 2:N_x-1
             if i == N_if1                  % interface 1
                 R(i,j) = ( ...
                    eps_r1*V(i+1,j) + ...
                    eps_r2*V(i-1,j) + ...
                    eps_av*V(i,j+1) + ...
                    eps_av*V(i,j-1))/(2*(eps_r1+eps_r2)) - V(i,j);
                 V(i,j) = V(i,j) + omega*R(i,j);
             elseif i == N_if2              % interface 2
                 R(i,j) = ( ...
                    eps_r2*V(i+1,j) + ...
                    eps_r1*V(i-1,j) + ...
                    eps_av*V(i,j+1) + ...
                    eps_av*V(i,j-1))/(2*(eps_r1+eps_r2)) - V(i,j);
                 V(i,j) = V(i,j) + omega*R(i,j);
             else                           % outside interfaces
                 R(i,j) = (V(i-1,j) + V(i+1,j) + V(i,j-1) + V(i,j+1))/4 - V(i,j);V(i,j) = V(i,j) + omega*R(i,j);
             end  
        end
    end
    Rmax = max(max(R));
    itnumber = itnumber + 1;
end

figure (1);
% plot the potencial countours
contour(V,50);
hold on;
% plot the vectors of electric field
[fx,fy] = gradient(-V);
quiver(fx,fy);
% visualize colors
colorbar
colormap winter
% add a horizontal line to visualize the planar interface 1
if1_up = yline(N_if1);
if1_up.Color = [0.4 0.4 0.4];
if1_up.LineStyle = '--';
if1_up.Label = ['\epsilon_1 = ' int2str(eps_r1) '\epsilon_0'];
if1_up.LabelVerticalAlignment = 'top';

if1_down = yline(N_if1);
if1_down.Color = [0.4 0.4 0.4];
if1_down.LineStyle = '--';
if1_down.Label = ['\epsilon_2 = ' int2str(eps_r2) '\epsilon_0'];
if1_down.LabelVerticalAlignment = 'bottom';

% add a horizontal line to visualize the planar interface2
if2_up = yline(N_if2);
if2_up.Color = [0.4 0.4 0.4];
if2_up.LineStyle = '--';
if2_up.Label = ['\epsilon_2 = ' int2str(eps_r2) '\epsilon_0'];
if2_up.LabelVerticalAlignment = 'top';

if2_down = yline(N_if2);
if2_down.Color = [0.4 0.4 0.4];
if2_down.LineStyle = '--';
if2_down.Label = ['\epsilon_1 = ' int2str(eps_r1) '\epsilon_0'];
if2_down.LabelVerticalAlignment = 'bottom';
hold off
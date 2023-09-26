%%
close all; clear; clc
%  -----------------------
% |                       |
% | --------------------- | <- IF2 = Ny/2+h/2
% | ////////eps_r//////// | 
% | --------------------- | <- IF1 = Ny/2-h/2
% |   ==w==   g   ==w==   |
%  -----------------------

%% Initialization
% universal constants
eps_0 = 8.85418782e-12;                         % permittivity of free space        
c = 299792458;                                  % speed of light
% mesh dimensions
Nx = 400;                                       % even number of nodes in x
Ny = 100;                                       % even number of nodes in y
% Ny/2 - h/2 - 1 > 20*h
h = min(2.5,round((Ny-2)/21));                  % dielectric height
if mod(h,2) == 1
    h = h + 1;
end
% Nx/2 - g/2 - w/2 > 20*h
w = min(58,round((Nx-40*h)/2.5));               % strips width
if w <= 0
    error('The ratio Nx/Ny is too low')
end
if mod(w,2) == 1
    w = w + 1;
end
g = round(2*w);                                 % gap between strips' centres
if mod(g,2) == 1
    g = g + 1;
end
IF1 = round(Ny/2-h/2);                          % 1st material interface
IF2 = round(Ny/2+h/2);                          % 2nd material interface
% electrostatic properties
V0 = 100;                                       % feeding potential
Vgnd = 0;                                       % ground potential
eps_r = 3.05;                                   % relative permittivity of the dielectric
er_avg = (1+eps_r)/2;
% initial values of potential across the mesh
V_init = zeros(Ny,Nx);                          % free nodes inside the structure
V_init(1,:) = Vgnd;                             % bottom edge of the bounding box
V_init(Ny,:) = Vgnd;                            % top edge of the bounding box
V_init(:,1) = Vgnd;                             % left edge of the bounding box
V_init(:,Nx) = Vgnd;                            % right edge of the bounding box
V_init(IF1-1,Nx/2-g/2-w/2:Nx/2-g/2+w/2) = V0;   % left strip (V0)
V_init(IF1-1,Nx/2+g/2-w/2:Nx/2+g/2+w/2) = -V0;  % right strip (-V0)
    
%% Calculation of potential -> charge -> capacity -> characteristic impedance
% Step 1: Obtain potential using the method of successive over-relaxation (SOR)
R = zeros(Ny,Nx);                               % improvements over the last iteration
% Step 2: Obtain charge using a discrete form of the Gauss theorem
% $\oint_{curve} \vec E \cdot \mathrm d \vec s = Q_{enc}$
% Positively oriented integration curve around the electrodes:
%  <-----D------
% | /////////// | ... IF1
% A  =w= g =w=  C
% |             |
%  ------B----->
k = 5;                                          % margin of the integration curve
left = Nx/2 - g/2 - w/2 - k;
bottom = IF1 - 1 - k;
right = Nx/2 - g/2 + w/2 + k;
top = IF1 - 1 + k;

tic
alpha_step = 0.01;
alpha_range = 1:alpha_step:1.99;                % examined values of alpha in [1,2)
alpha_range = 1.94;                             % single value for debugging with nice fig2 (1.7)
alpha_iter = ones(1,length(alpha_range));       % number of iterations for given alpha
for alpha = alpha_range                         % calculate everything for different alphas
    Rmax = V0;                                  % largest point improvement across the mesh
    iter = 0;                                   % number of iterations for current alpha
    V = V_init;                                 % initialize potential values
    while Rmax > 1e-3                           % relaxation until improvement is insignificant
        iter = iter + 1;
        % Step 1: SOR method
        for j = 2:Nx-1
            for i = 2:Ny-1
                % nodes on interface 1
                if i == IF1
                    R(i,j) = (eps_r*V(i+1,j)+1*V(i-1,j)+er_avg*V(i,j-1)+er_avg*V(i,j+1))/4/er_avg-V(i,j); 
                % nodes on interface 2
                elseif i == IF2
                    R(i,j) = (1*V(i+1,j)+eps_r*V(i-1,j)+er_avg*V(i,j-1)+er_avg*V(i,j+1))/4/er_avg-V(i,j);
                % nodes in strips
                elseif i == IF1-1 && ((j >= Nx/2-g/2-w/2 && j <= Nx/2-g/2+w/2) || (j >= Nx/2+g/2-w/2 && j <= Nx/2+g/2+w/2))
                    R(i,j)= R(i,j);
                    V(i,j)= V(i,j);
                % free nodes
                else
                    R(i,j) = (V(i+1,j)+V(i-1,j)+V(i,j-1)+V(i,j+1))/4-V(i,j);
                end
                V(i,j) = V(i,j)+alpha*R(i,j);
            end
        end
        Rmax = max(max(R));
        % Step 2: numerical evaluation of the integral along the curve
        qA = 0; qB = 0; qC = 0; qD = 0;
        for i=top:1:bottom
            if i < IF1
                qA = qA + eps_r*(V(i,left+1) - V(i,left-1))/2;
                qC = qC + eps_r*(V(i,right-1) - V(i,right+1))/2;
            elseif i == IF1
                qA = qA + er_avg*(V(i,left+1) - V(i,left-1))/2;
                qC = qC + er_avg*(V(i,right-1) - V(i,right+1))/2;
            else
                qA = qA + (V(i,left+1) - V(i,left-1))/2;
                qC = qC + (V(i,right-1) - V(i,right+1))/2;
            end
        end
        for j=left:1:right
                qB = qB + (V(bottom+1,j) - V(bottom-1,j))/2;
                qD = qD + eps_r*(V(top-1,j) - V(top+1,j))/2;
        end
        Q(iter) = eps_0*(qA + qB + qC + qD);    % charge per unit length on both electrodes
    end
    C = Q/(V0-Vgnd);                            % capacitance per unit length of the structure
    Z = 1./(c*C);                               % charactersitic impedance Z0
    alpha_iter(alpha_range==alpha) = iter;      % save number of iterations for current alpha
end
toc

%% Results
% Resultant electric field distribution within task geometry
% Important note: Results are saved for the last value of alpha since all
% results should be representative
figure (1);
% plot the potential countours
contour(V,20);
hold on;
% plot the vectors of electric field
[Ex,Ey] = gradient(-V);
quiver(Ex,Ey);
% plot settings
xticklabels(round((xticks-Nx/2)./h))
yticklabels(round((yticks-Ny/2)./h))
xlabel('x [mm]','FontSize',18)
ylabel('y [mm]','FontSize',18)
colorbar
% colormap winter
% colormap jet
title(['Z = ' num2str(round(Z(end),1)) ...
    ' \Omega, C = ' num2str(round(C(end)*1e9,3)) ' nF/m'],'FontSize',18)
subtitle(['w = ' num2str(round(w/h,2)) ' mm, g = ' num2str(round(g/h,2)) ' mm'],'FontSize',18)
% interface 1
if1_up = yline(IF1);
if1_up.Color = [0.4 0.4 0.4];
if1_up.LineStyle = '--';
% if1_up.Label = ['\epsilon_2 = ' num2str(eps_r) '\epsilon_0'];
if1_up.LabelVerticalAlignment = 'top';
if1_up.LabelHorizontalAlignment = 'left';
if1_down = yline(IF1);
if1_down.Color = [0.4 0.4 0.4];
if1_down.LineStyle = '--';
% if1_down.Label = '\epsilon_1 = \epsilon_0';
if1_down.LabelVerticalAlignment = 'bottom';
if1_down.LabelHorizontalAlignment = 'left';
% interface 2
if2_up = yline(IF2);
if2_up.Color = [0.4 0.4 0.4];
if2_up.LineStyle = '--';
% if2_up.Label = '\epsilon_1 = \epsilon_0';
if2_up.LabelVerticalAlignment = 'top';
if2_up.LabelHorizontalAlignment = 'left';
if2_down = yline(IF2);
if2_down.Color = [0.4 0.4 0.4];
if2_down.LineStyle = '--';
% if2_down.Label = ['\epsilon_2 = ' num2str(eps_r) '\epsilon_0'];
if2_down.LabelVerticalAlignment = 'bottom';
if2_down.LabelHorizontalAlignment = 'left';
hold off

% Evolution of the characteristic impedance with the number of iterations
figure(2)
plot(Z)
axis tight
ylabel('Z_0 [\Omega]','FontSize',18)
xlabel('Number of iterations [-]','FontSize',18)
ylim([0 100])
% Z0 = 50 line
if1_up = yline(50);
if1_up.Color = [0.4 0.4 0.4];
if1_up.LineStyle = '--';

%
figure(3)
plot(Ey(IF1-1,Nx/2-g/2-w/2:Nx/2-g/2+w/2))
axis tight
xticklabels(round((xticks-w/2)./h))
xlabel('Strip-centered x-coordinate [mm]','FontSize',18)
ylabel('Charge density','FontSize',18)

% Compare numerically optimal value of alpha with the analytical solution
r = cos(pi/Nx) + cos(pi/Ny);
alpha_optimal = (8-sqrt(64-16*r^2))/r^2;
disp(['Analytical solution for the optimal value of alpha: ' num2str(round(alpha_optimal,2))])
disp(['Numerically determined optimal value of alpha: ' num2str(alpha_range(alpha_iter==min(alpha_iter)))])

% Plot the dependecy of the number of iterations for given alpha 
figure(4)
plot(alpha_range,alpha_iter)
axis tight
ylabel('Number of iterations [-]','FontSize',18)
xlabel('SOR coefficient \alpha [-]','FontSize',18)

% Full run (alpha = 1:0.01:1.9) old results:
% Elapsed time is 3969.8 seconds.
% Analytical solution for the optimal value of alpha: 1.96
% Numerically determined optimal value of alpha: 1.91

% Full run (alpha = 1:0.01:1.9) new results:
% Elapsed time is 70.186527 seconds.
% Analytical solution for the optimal value of alpha: 1.96
% Numerically determined optimal value of alpha: 1.94

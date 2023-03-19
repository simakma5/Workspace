%%
close all; clear; clc

%% Calcualtion of charactersitic impedance of a square coaxial cable with a strip inner conductor
% using Successive over Relaxation (SOR) method

c = 3*1e8; % speed of light
eps = 8.845*1e-12; % permitivity of vacuum

%the size of the computational domain
Nx=50; Ny=Nx; %number of nodes, must be even number!

% set potencials at fixed nodes (electrodes)
V0=100; Vgnd=0; 
 
Nd1=round(1*Ny/4); %row number of #1 dielectric interface
Nd2=round(3*Ny/8); %row number of #1 dielectric interface
er1 = 1;
er2 = 1;
eavr = (er1+er2)/2;
 
% set initial potencial guess at free (inner) nodes equal to the average of fixed values
for i=2:Ny-1
    for j=2:Nx-1
        V(i,j)=(V0+Vgnd)/2;
    end
end

% set potencials at fixed nodes (top and bottom sides)
V(1,:) = Vgnd;
V(Ny,:) = Vgnd; %Note. Matlab displays the lowest row on the top of the picture.    

% set potencials at fixed nodes (left and right sides)
V(:,1) = Vgnd;
V(:,Nx) = Vgnd;

% inner point electrode
% V((Ny-1)/2+1,(Nx-1)/2+1) = V0;

% inner strip electrode
W = 20; % strip width, even number!
V(Ny/2,Nx/2-W/2:Nx/2+W/2) = V0;
    
%% update potencial values at free nodes by Successive Over Relaxation (SOR)
omega=1.7; %omega - SOR coefficient is choosen from the interval <1;2)
tres = 0.001;
Rmax = V0;
niter = 0;  %number of iterations
while (Rmax > tres)
    for i = 2:Ny-1
        for j = 2:Nx-1
             if (i ~= Nd1) && (i ~= Nd2) %nodes out of the interface
                R(i,j)= 1/4*(V(i+1,j)+V(i-1,j)+V(i,j-1)+V(i,j+1))-V(i,j);
                V(i,j)= V(i,j)+omega*R(i,j);
                V(Ny/2,Nx/2-W/2:Nx/2+W/2) = V0;
             elseif i == Nd1 %nodes on the interface
                R(i,j)= 1/4*(er2*V(i+1,j)+er1*V(i-1,j)+eavr*V(i,j-1)+eavr*V(i,j+1))/eavr-V(i,j); 
                V(i,j)= V(i,j)+omega*R(i,j);
                V(Ny/2,Nx/2-W/2:Nx/2+W/2) = V0;
             elseif i == Nd2 %nodes on the interface
                R(i,j)= 1/4*(er1*V(i+1,j)+er2*V(i-1,j)+eavr*V(i,j-1)+eavr*V(i,j+1))/eavr-V(i,j); 
                V(i,j)= V(i,j)+omega*R(i,j);
                V(Ny/2,Nx/2-W/2:Nx/2+W/2) = V0;
             end  
        end
    end
    Rmax = max(max(R));
    niter = niter + 1;
end

%% charge, capacity, characteristic impedance
k = 2; % distance of the integration curve from the strip, try to set higher k = 3,4,5, ..

left = Nx/2 - W/2 - k;
right = Nx/2 + W/2 + k;
up = Ny/2 - k;
down = Ny/2 + k;

QQ = 0;

% the cross section of the coaxial cable, 
% A,B,C,D denotes the parts of the integration curve
% 
%  -----B-----
% |           |
% A   --W--   C
% |           |
%  -----D-----

% tic;
% % vertical parts of the integration curve A and C
% for i=up:1:down % e.g. for Ny=50, k=2: i = 23:1:27
%          QQ = QQ + (V(i,left+1) - V(i,left-1))/2; %A 
%          QQ = QQ + (V(i,right-1) - V(i,right+1))/2; %C
% end
% 
% % horizontal parts of the integration curve B and D
% for j=left:1:right
%          QQ = QQ + (V(up+1,j) - V(up-1,j))/2; %D, 24-22
%          QQ = QQ + (V(down-1,j) - V(down+1,j))/2; %B, 26-28
% end
% toc;

%% Vectorized calculation. 
% Check computational time consumed if "for" cycle and vectorized versions are used for different size of the computational domain 

tic;
% vertical parts of the integration curve A and C
Ql = V(up:down,left+1) - V(up:down,left-1);
Qr = V(up:down,right-1) - V(up:down,right+1); 
QV = sum(Ql)/2+sum(Qr)/2;

% horizontal parts of the integration curve B and D
Qt = V(up+1,left:right) - V(up-1,left:right);
Qb = V(down-1,left:right) - V(down+1,left:right); 
QH = sum(Qt)/2+sum(Qb)/2;
  
QQ = QV + QH;
toc;

Q = eps*QQ; % charge
C = Q/(V0-Vgnd);   % capacity
Z0 = 1/(c*C); % charactersitic impedance Z0

% Z0 analytically, valid only for vacuum (air dielectric)
% according to Bohumil Šimíček, Antény pro televizní a rozhlasové vysílání na VKV, Nadas, Praha, 1989, str. 33.
Z0an = 130*log10(2.16*Nx/W); % valid for Nx/W > 2

%% figures
figure (1);
contour(V,10);  %plot the potencial countours
hold on;
[fx,fy] = gradient(-V);
quiver(fx,fy);	%plot the vectors of electric field
colorbar;

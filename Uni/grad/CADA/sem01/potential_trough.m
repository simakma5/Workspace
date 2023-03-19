% set potencials at fixed nodes (electrodes)
V0=100; Vgnd=0; 
% number of nodes 
Nx=10; Ny=Nx; 

%set initial guess of free (inner) nodes equal to average of fixed values
V = zeros(Ny,Nx);
for i = 2:Ny-1
    for j = 2:Nx-1
        V(i,j)=(V0+Vgnd)/2;
    end
end
%set potencials at fixed nodes (top side), V=V0
V(Ny,:) = V0;
% set potencials at fixed nodes (bottom, left and right sides), V=Vgnd
V(1,:) = Vgnd;
V(:,1) = Vgnd;
V(:,Nx) = Vgnd;
% update potencial values at free nodes by Successive Over Relaxation (SOR)
omega=1.5; % omega - SOR coefficient is choosen from the interval <1;2)
tresh = V0/10e6;
Rmax = V0;
iter = 0;
R = zeros(Ny,Nx);
while Rmax > tresh
    for i = 2:Ny-1
       for j = 2:Nx-1
           R(i,j) = (V(i-1,j) + V(i+1,j) + V(i,j-1) + V(i,j+1))/4 - V(i,j);
           V(i,j) = V(i,j) + omega*R(i,j);
       end
    end
Rmax = max(max(R));
iter = iter + 1;
end

%plot of electric field distribution
figure (1)
contour(V,20);
hold on;
[fx,fy] = gradient(-V);
quiver(fx,fy);		
colorbar;
disp(iter);

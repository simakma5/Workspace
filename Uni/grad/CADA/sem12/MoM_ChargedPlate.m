close all; clear; clc

%% Calculation of the charge density at the plate using the Method of Moments (MoM)
% square number of cells (= number of unknowns, i.e. basis functions)
N = 1e4;
% plate size 2a x 2a
a = 1;
% cell size 2b x 2b
b = a/sqrt(N);

eps = 8.85e-12;

% Assignment of cell centers into x, y vectors of the size (1xN)
x = zeros(1,N);
y = zeros(1,N);
for i = 1:sqrt(N)
   for j = 1:sqrt(N)  
       x((i-1)*sqrt(N)+j) = a-b-(j-1)*2*b;
       y((i-1)*sqrt(N)+j) = a-b-(i-1)*2*b;
   end
end

% Elements of the Amn matrix
A = zeros(N,N);
for m = 1:N
   for n = 1:N
       if m == n
          A(m,n) = 0.8814*2*b/(pi*eps); 
       else
          A(m,n) = b^2/(pi*eps*sqrt((x(m)-x(n))^2+(y(m)-y(n))^2));
       end
   end
end
alpha = A\ones(N,1);

% Charge density at the plate
chargeDensity = zeros(sqrt(N),sqrt(N));
for i = 1:sqrt(N)
    for j = 1:sqrt(N)
        chargeDensity(i,j) = alpha((i-1)*sqrt(N)+j);           
    end
end

% Plot results
figure(2)
mesh(chargeDensity)
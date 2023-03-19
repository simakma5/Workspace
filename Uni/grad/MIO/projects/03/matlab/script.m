%% Ideal Directional Coupler
Zv = 50;
C = 3;

Zve = Zv*sqrt((1+10^(-C/20))/(1-10^(-C/20)));
Zvo = Zv*sqrt((1-10^(-C/20))/(1+10^(-C/20)));

% assert Zcheck == Zv
Zcheck = sqrt(Zve*Zvo);
Zdiff = round(abs(Zve-Zvo),2);

%% Coupled MSLines length
% from TXLine
lambda_even = 17.2;
lambda_odd = 18.6;

lambda = 2*lambda_odd*lambda_even/(lambda_odd+lambda_even);
L = lambda/4;

%% Two Branch Coupler
Zv = 50;
% C = 3; implies relations in lectures

Zv1 = Zv;
Zv2 = Zv/sqrt(2);


%% Three Branch Coupler
Zv = 50;
% C = 3; implies relations in lectures

Zv1 = Zv/(sqrt(2)-1);
Zv2 = Zv/sqrt(2);
Zv3 = Zv2;

%% Ideal Rat Race Coupler
Zv = 50;
% C = 3; implies relations in lectures

Zv1 = sqrt(2)*Zv;
Zv2 = Zv1;

%% Uncompensated Wilkinson Power Divider
Zv = 50;
% C = 3; implies relations in lectures

Zv1 = sqrt(2)*Zv;
Zv2 = Zv1;

%% Compensated Wilkinson Power Divider
Zv = 50;
% C = 3; implies relations in lectures

Zv1 = 2^(-1/4)*Zv;
Zv2 = 2^(+1/4)*Zv;
R = 2*Zv;

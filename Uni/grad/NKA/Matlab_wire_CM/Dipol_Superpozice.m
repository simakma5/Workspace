clear all; close all; clc;

%%% Vypocet vstupni impedance dipolu

mi        = 4*pi*1e-7;    % permeability of vacuum
epsilon   = 8.85e-12;     % permittivity of vacuum
c         = 3e8;

N         = 101;          % pocet segmentu

L         = 0.1;          % celkova delka dipolu [m]
a         = 0.1e-3;       % polomer dipolu [m]

feed      = floor(N/2+1); % napajeci segment (stred dipolu), muze byt obecny

Ei        = zeros(N,1);   % budici vektor [Ei]

%Ei(feed)  = 10;           % el. pole na stredovem napajecim segmentu
Ei        = ones(N,1);   % rovinna vlna
%Ei(25)    = 10;           % napajeni v 1/4 delky dipolu
%Ei(10) = 5;     % ruzne moznosti buzeni




 % priklad zobrazeni proudu na dane frekvenci
 f = 5.0e9;

 [Z,delta] = Zmatrix(f,L,N,a,0);
 Jmom      = inv(Z)*Ei;

 
 [We,Wm,Pr] = Wstored11(Jmom,N,delta,f,a) 
  
 % priklad rozkladu na charakteristicke mody
 R = real(Z);
 X = imag(Z);
 [V,D] = eig(X,R,'qz');
 D = diag(D);
 
 
 % Superpozice
 
 midxstart = 12 % index od ktereho zacinaji mody
 modu      = 4 % kolik modu secist (min 4)

 m = 1;
 
 for p = midxstart:-1:midxstart-modu+1;
 J(:,m)  = V(:,p)./sqrt(V(:,p)'*R*V(:,p));
 MA(m)   = 1/(1+i*D(p));
 ME(m)   = sum(Ei.*J(:,m));
 A(m)    = MA(m)*ME(m);
 Js(:,m) = A(m).*J(:,m);
 m=m+1;
 end;

figure;
plot(J(:,1),'r','LineWidth',2);
hold on;
plot(J(:,2),'g','LineWidth',2);
hold on;
plot(J(:,3),'b','LineWidth',2);
hold on;
plot(J(:,4),'k','LineWidth',2);
title('Prvni 4 normovane charakteristicke proudy');
xlabel('Segment [-]');
ylabel('I [A]');
legend('J1','J2','J3','J4');
grid on;


 
 Jtotal = sum(Js(:,:)')';
 
 figure;
 subplot(2,1,1);
 plot(1e3.*real(Jmom),'b','LineWidth',2);
 hold on;
 plot(1e3.*imag(Jmom),'r','LineWidth',2);
 hold on;
 plot(1e3.*abs(Jmom),'k','LineWidth',2);
 legend('Re(I)','Im(I)','|I|');
 xlabel('segment [-]');
 ylabel('I [mA]');
 title('Proud ziskany primou inverzi [Z]*[Ei]');
 
 subplot(2,1,2); 
 plot(1e3.*real(Jtotal),'b','LineWidth',2);
 hold on;
 plot(1e3.*imag(Jtotal),'r','LineWidth',2);
 hold on;
 plot(1e3.*abs(Jtotal),'k','LineWidth',2);
 legend('Re(I)','Im(I)','|I|');
 xlabel('segment [-]');
 ylabel('I [mA]');
 title('Proud ziskany superpozici char. modu');

 
 figure;
 subplot(3,1,1);
 stem(abs(MA),'Color','red','LineWidth',2);
 title('Modalni amplituda |MA|');
 subplot(3,1,2);
 stem(ME,'Color','blue','LineWidth',2);
 title('Modalni excitacni koeficient ME');
 subplot(3,1,3);
 stem(abs(A),'Color','black','LineWidth',2);
 title('Expanzni koeficienty pro celkovy proud |MA*ME|');
MA
ME
A

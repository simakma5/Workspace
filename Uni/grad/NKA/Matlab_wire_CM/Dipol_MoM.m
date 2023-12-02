clear all; close all; clc;

%%% Vypocet vstupni impedance dipolu a cinitelu jakosti
% update 10.11.2015

mi        = 4*pi*1e-7;    % permeability of vacuum
epsilon   = 8.85e-12;     % permittivity of vacuum
c         = 3e8;

N         = 61;          % pocet segmentu

L         = 0.1;          % celkova delka dipolu [m]
a         = 0.1e-3;       % polomer dipolu [m]

%feed      = floor(N/2+1); % napajeci segment (stred dipolu), muze byt obecny

feed = 25  % napajeni v 1/4 dipolu

Ei        = zeros(N,1);   % budici vektor [Ei]
Ei(feed)  = floor(N/2);           % el. pole na napajecim segmentu


% frekvencni osa
fstart = 1.0e9; fend = 5e9; fstep = 50e6;
fosa = (fstart:fstep:fend)./1e9;

m = 1;
df = 1e6;
for f = fstart:fstep:fend;

  [Z,delta]           = Zmatrix(f,L,N,a,0); % vypocet impedancni matice
  J                   = inv(Z)*Ei;         % vypocet vektoru proudu [J] = [Z]^{-1}[Ei]
  Jatfeed(m)          = J(feed);           % proud na napajecim segmentu pro vypocet Zin
  Zin(m)              = Ei(feed)/J(feed); % vypocet vstupni impedance
  [We(m),Wm(m),Pr(m)] = Wstored11(J,N,delta,f,a);
  Qmax(m)             = 2*2*pi*f*max(We(m),Wm(m))/Pr(m);
  Q(m)                = 2*pi*f*(We(m)+Wm(m))/Pr(m);

  
  % Qz, Qx
  [Zplus,delta]       = Zmatrix(f+df,L,N,a,0);
  Jplus               = inv(Zplus)*Ei;  
  [Zminus,delta]      = Zmatrix(f-df,L,N,a,0);
  Jminus              = inv(Zminus)*Ei;  
  Zinplus             = Ei(feed)/Jplus(feed);
  Zinminus            = Ei(feed)/Jminus(feed);
  R0                  = real(Zin(m));
  dZ                  = Zinplus-Zinminus;
  dX                  = imag(Zinplus-Zinminus);
  dom                 = 2*pi*2*df;
  Qz(m)               = (2*pi*f)*abs(dZ/dom)/(2*R0);
  Qx(m)               = (2*pi*f)*(dX/dom)/(2*R0);
  
  m = m +1;
 end;
 
 % Return Loss
 Z0 = 73;
 RL = 20*log10(abs((Zin-Z0)./(Zin+Z0)));

 % vykresleni Zin(f)
 figure;
 hold on;
 plot(fosa,real(Zin),'k','LineWidth',2);
 plot(fosa,imag(Zin),'r','LineWidth',2);
 xlabel('f [GHz]');
 ylabel('Z [\Omega]');
 legend('Rin','Xin');
 grid on;
 
  figure;
 hold on;
 plot(fosa,Q,'LineWidth',2);
 plot(fosa,Qmax,'LineWidth',2);
 plot(fosa,Qz,'r');
 plot(fosa,Qx,'k--');
 legend('Q','Q (max[W_m,W_e]','Q_Z','Q_X');
 xlabel('f [GHz]');
 ylabel('Q [-]');
 grid on;
 
 figure;
 plot(fosa,RL,'LineWidth',2);
 xlabel('f [GHz]');
 ylabel('RL [dB]');
 grid on;
 
 % priklad zobrazeni proudu na dane frekvenci
 f = 4e9;
 
 [Z,delta] = Zmatrix(f,L,N,a,0);
 J         = inv(Z)*Ei;
  
 figure;
 hold on;
 plot(real(J),'b','LineWidth',2);
 plot(imag(J),'r','LineWidth',2);
 legend('Re(I)','Im(I)');
 xlabel('segment [-]');
 ylabel('I [A]');
 title('Celkovy proud na 4GHz')
 grid on;
 
 
 % priklad rozkladu na charakteristicke mody
 R = real(Z);
 X = imag(Z);
 [V,D] = eig(X,R);
 D = diag(D);
 
 % pozn. na teto frekvenci jsou indexy prvnich modu 7,6,5,4
 
 idx1  = 7;
 idx2  = 6;
 idx3  = 5;
 
 mode1 = V(:,idx1);
 mode2 = V(:,idx2);
 mode3 = V(:,idx3);
 
 % normovani na Pr=1W
 mode1N = sqrt(2).*V(:,idx1)./sqrt(V(:,idx1)'*R*V(:,idx1));
 mode2N = sqrt(2).*V(:,idx2)./sqrt(V(:,idx2)'*R*V(:,idx2));
 mode3N = sqrt(2).*V(:,idx3)./sqrt(V(:,idx3)'*R*V(:,idx3));
 
 % nenormovane proudy
 figure;
 hold on;
 plot(mode1,'r','LineWidth',2);
 plot(mode2,'g','LineWidth',2);
 plot(mode3,'b','LineWidth',2);
 legend('J1','J2','J3');
 xlabel('segment [-]');
 ylabel('I [A]');
 title('Nenormovane char. proudy po vystupu z eig');
 grid on;
 
 
 % normovane na Pr=1W
 figure;
 hold on;
 plot(mode1N,'r','LineWidth',2);
 plot(mode2N,'g','LineWidth',2);
 plot(mode3N,'b','LineWidth',2);
 legend('J1','J2','J3');
 xlabel('segment [-]');
 ylabel('I [A]');
 title('Normovane char. proudy Pr=1W po vystupu z eig')
 grid on;
 
 disp('Q jednotlivych char. proudu'); 
 [We1,Wm1,Pr1] = Wstored11(mode1,N,delta,f,a);
 [We2,Wm2,Pr2] = Wstored11(mode2,N,delta,f,a);
 [We3,Wm3,Pr3] = Wstored11(mode3,N,delta,f,a);
 Q1 = 2*2*pi*f*max(We1,Wm1)/Pr1
 Q2 = 2*2*pi*f*max(We2,Wm2)/Pr2
 Q3 = 2*2*pi*f*max(We3,Wm3)/Pr3





 
 
 

  
  

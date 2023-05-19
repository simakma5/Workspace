%Code: dipole.m		Courtesy: S. Vivek
% current distribution and input impedance of a center fed wire dipole using
% Poklington's integral equation is solved.
% pulse expansion and point matching employed.
% provision for magnetic-frill source model
% symmetry about the middle is used to reduce the size.
% the typical input data: lambda = 1, radius = 0.005, L = 0.47

close all;
clc;
clf;

syms z
c=3e8; 
ep0=8.85e-12;  
lambda = input('enter wavelength in meter: ');
l0 = lambda;
f=c/l0;   	    % frequency
w=2*pi*f;   	% angular frequency
Vo=1;           % excitation voltage
radius = input('enter the radius of dipole in wavelength: ');
a=radius*l0;
L = input('enter the length of dipole in wavelength: ');
l = L*l0;
k=2*pi/l0;     	% wave number
N1=[11 21 31 51 71 81];		% number of basis functions or segments, odd number only
Ii=zeros(1,length(N1)); 
za = Ii;
for x=1:length(N1)
    N=N1(x);
    h=l/N;     			% segment length
    K = zeros(N);
    zm=((1:N)-0.5).*(l/N)-(l/2);  
    Is=(N+1)/2;
    V=zeros(N,1);
    F=inline((exp(-1i*k*(a^2+(zm(1)-z).^2).^0.5)./(4*pi*(a^2+(zm(1)-z).^2).^2.5))...
                .*(((1+1i*k*(a^2+(zm(1)-z).^2).^0.5).*(2*(a^2+(zm(1)-z).^2)-3*a^2))...
                +((k*a)^2*(a^2+(zm(1)-z).^2))));       % Poklington's Equation

        for jn=1:N
            K(1,jn)=quad(F,zm(jn)-(h/2),zm(jn)+(h/2));     % Numerical integration 
        end
        K(1,1) = K(1,1)/2; 
        for in = 2:N            	% use of Symmetry
            for jn = in:N
                K(in,jn) = K(in-1,jn-1);
            end
        end
        K = K+transpose(K);       
        
        r1=sqrt(zm.^2+(a)^2);
        r2=sqrt(zm.^2+(2.3*a)^2);
	%Magnetic Frill generator
%      V=transpose(-i*w*ep0*Vo/2/log(2.3)*(exp(-i*k*r1)./r1-exp(-i*k*r2)./r2));       
 	%Delta gap generator
      V(Is) = -1i*w*ep0*Vo/h; 
        
    alpha= K\V; 	% current distribution
    Ii(1:N,x)=alpha;
    za(1:N,x)=zm;
    Zin(x)=1/alpha(Is);     % input impedance Zin=Vo/Ii
end

figure(1);
plot(za(1:N1(1),1),abs(Ii(1:N1(1),1)),'r-',za(1:N1(2),2),abs(Ii(1:N1(2),2)),'b-',za(1:N1(3),3),...
    abs(Ii(1:N1(3),3)),'k-',za(1:N1(4),4),abs(Ii(1:N1(4),4)),'k--', za(1:N1(5),5),abs(Ii(1:N1(5),5)),'r--', za(1:N1(6),6),abs(Ii(1:N1(6),6)),'g--');
legend('N=11','N=21','N=31','N=51','N=71','N=81');
title('current distribution on a center fed dipole for various value of N');
xlabel('z/\lambda');
ylabel('Current distribution, abs(I)');

figure(2)
subplot(2,1,1),plot(N1,real(Zin),'k-*');grid;
title('input impedance of the dipole for various value of N');
xlabel('number of basis functions, N ');
ylabel('Rin');
subplot(2,1,2),plot(N1,imag(Zin),'k-*'); grid;
xlabel('number of basis functions, N ');
ylabel('Xin');
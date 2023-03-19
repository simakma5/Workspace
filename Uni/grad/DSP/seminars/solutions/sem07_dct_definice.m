% Cviceni na Diskretni Kosinovou Transformaci
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRVNI CAST - definice DCT a souvislost DFT a DCT

clear all;
% close all ;


x=[10 7 4 1 -2 -5 -8 -5 -2 1] ;
% x=1:10;
N=length(x);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DCT-1 - 2N-2 extension (Oppenheim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)

subplot(221)
n=0:(N-1) ;
stem(n,x) ;
title('Original sequence')

% % DCT-1 podle Oppenheima (2N-2 extension) 
XC1 = dctxc1(x) ;

kk=0:N-1;
subplot(222); h=plot(kk,XC1,'b.-','MarkerSize',20); title('DCT-1 (Oppenheim) - dctxc1.m')

x01=idctxc1(XC1);
subplot(221)
hold on; stem(n,x01,'r') ; hold off

xx=[x fliplr(x(2:N-1))];   % (2N-2)-point even extension

subplot(223);
nn=0:(2*N-3) ;
stem(nn,xx) ;
title('2N-2 extension')

kkk=0:2*N-3;
XX=real(fft(xx)) ;   % real is used to remove a parazit imaginary component appearing due to finite precision computation
subplot(224); plot(kkk,XX,'b.-','MarkerSize',20); hold on ; plot(kk,XC1,'r.-','MarkerSize',20); hold off
title('DFT - 2N-2 even extension')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  DCT-2 - 2N extension (Oppenheim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


beta=[1/2 ones(1,N-1)]';


figure(2)
subplot(221)
n=0:(N-1) ;
stem(n,x) ;
title('Original sequence')

XC2=dctxc2(x);

kk=0:N-1;
subplot(222);  plot(kk,XC2,'b.-','MarkerSize',20); 
title('DCT-2 (2N extension) - dctxc2.m');

x02=idctxc2(XC2);
subplot(221)
hold on; stem(n,x02,'r') ; hold off


xx=[x fliplr(x(1:N))];   % 2N-point even extension (dct in MATLAB)


subplot(223);
nn=0:(2*N-1) ;
stem(nn,xx) ;
title('2N extension')

XX=fft(xx); 
k=1:2*N;
W=exp(-j*pi*(k-1)/(2*N));
XXX=real(XX.*W); 

kkk=0:2*N-1;
subplot(224)
plot(kkk,XXX,'b.-','MarkerSize',20)
hold on
plot(kk,XC2,'r.-','MarkerSize',20)
hold off

% x022=idctxc2(XXX)
% subplot(221)
% hold on; stem(n,x022,'g') ; hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  DCT-2 - 2N extension (MATLAB - dct function)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(3)
subplot(221)
n=0:(N-1) ;
stem(n,x) ;
title('Original sequence')

% DCT in MATLAB (difference in scaling coeff -> multipl. by sqrt(N))
% XC2=sqrt(N*2)*dct(x);
XDCT=dct(x);

kk=0:N-1;
subplot(222);  plot(kk, XDCT,'b.-','MarkerSize',20); 
title('Matlab DCT (2N extension) - function dct()');


x02=idct(XDCT);
subplot(221)
hold on; stem(n,x02,'r') ; hold off


xx=[x fliplr(x(1:N))];   % 2N-point even extension (dct in MATLAB)


subplot(223);
nn=0:(2*N-1) ;
stem(nn,xx) ;
title('2N extension')

XX=fft(xx); 
k=1:2*N;
W=exp(-j*pi*(k-1)/(2*N));
XXX=real(XX.*W); 

kkk=0:2*N-1;
subplot(224)
plot(kkk,XXX,'b.-','MarkerSize',20)
hold on
plot(kk,XDCT,'r.-','MarkerSize',20)
hold off
ylim([-max(XXX) max(XXX)])

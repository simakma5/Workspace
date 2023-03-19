% Demo KLT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
close all; clear all

% Setup of directory with signals for different OS
if exist('/home/pollak/vyuka/DSP/m'),  % it is the Linux OS
    addpath /home/pollak/vyuka/DSP/m
end
if exist('/home/pollak/vyuka/DSP/signaly'),  % it is the Linux OS
    addpath /home/pollak/vyuka/DSP/signaly
end

if exist('H:\VYUKA\DSP\m'),
    addpath H:\VYUKA\DSP\m
end
if exist('H:\VYUKA\DSP\signaly'),
    addpath H:\VYUKA\DSP\signaly
end

if exist('C:\usr\pollak\vyuka\DSP\m')
    addpath C:\usr\pollak\vyuka\DSP\m
end
if exist('C:\usr\pollak\vyuka\DSP\signaly')
    addpath C:\usr\pollak\vyuka\DSP\signaly
end

% Ortogonality of eigen vectors

X=[1 0.5; 0.5 1];

[V, D] = eig(X); 

disp('Original matrix')
disp(X)
disp('Eigenvectors')
disp(V)
disp('Eigenvvalues in diagonal matrix')
disp(D)


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Karhunen-Loeve Transform
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% signal for the analysis


% one sinusoidal
fs=8000;
f=821 ;
s1=sin(2*pi*f*(1:2000)./fs) ;

% two sinusoidals
fs=8000;
f1=532 ;
f2=940 ;
s2=sin(2*pi*f1*(1:2000)./fs)+sin(2*pi*f2*(1:2000)./fs) ;

% voiced speech frame
s3=loadbin('vm0.bin');

s1=s1(:); 
s2=s2(:); 
s3=s3(:); 
slen1=length(s1);
slen2=length(s2);
slen3=length(s3);

% Addition of additive Gaussian noise
% scale = 0 ;
scale = 0.001 ;
% scale = 0.1 ;
s1=s1+scale*randn(slen1,1) ;
s2=s2+scale*randn(slen2,1) ;
s3=s3+scale*randn(slen3,1) ;


% Computation of covariance matrix from 20 realizations

wnum=20;
wlen=200;
wstep=50;

for i =1:wnum,
    SS1(i,:)=s1((i-1)*wstep+1:(i-1)*wstep+wlen) ;
    SS2(i,:)=s2((i-1)*wstep+1:(i-1)*wstep+wlen) ;
    SS3(i,:)=s3((i-1)*wstep+1:(i-1)*wstep+wlen) ;
end

RSS1=SS1'*SS1;
RSS2=SS2'*SS2;
RSS3=SS3'*SS3;

[V1, D1] = eig(RSS1); 
[V2, D2] = eig(RSS2); 
[V3, D3] = eig(RSS3); 

% Check of eigenvector ortogonality
disp(' ')
disp(['Scalar product of V3(:,10) & V3(:,20) = ' num2str(V3(:,10)'*V3(:,20))])
disp(['Scalar product of V3(:,120) & V3(:,140) = ' num2str(V3(:,120)'*V3(:,140))])
disp(['Scalar product of V3(:,150) & V3(:,150) = ' num2str(V3(:,150)'*V3(:,150))])

% Creation of transformation matrix of KLT

WKLT1=flipud(V1');
WKLT2=flipud(V2');
WKLT3=flipud(V3');

figure(1)
subplot(311)
plot(s1)
title('SIGNAL s1 - whole')

subplot(312)
plot(SS1(7,:)) ;
title('SIGNAL s1 - whole')

subplot(313)
stem(diag(D1)) ;
title('SIGNAL s1 - eigenvalues of cov.matrix')

figure(2)
subplot(411)
plot(WKLT1(1,:))
title('SIGNAL s1 - 1st basis (eigen)vector (most significant)')

subplot(412)
plot(WKLT1(2,:)) ;
title('2nd basis (eigen)vector')

subplot(413)
plot(WKLT1(3,:)) ;
title('3rd basis (eigen)vector')

subplot(414)
plot(WKLT1(4,:)) ;
title('4th basis (eigen)vector')

% %%%%%%%%%%%%

figure(3)
subplot(311)
plot(s2)
title('SIGNAL s2 - whole')

subplot(312)
plot(SS2(7,:)) ;
title('SIGNAL s2 - whole')

subplot(313)
stem(diag(D2)) ;
title('SIGNAL s2 - eigenvalues of cov.matrix')

figure(4)
subplot(411)
plot(WKLT2(1,:))
title('SIGNAL s2 - 1st basis (eigen)vector (most significant)')

subplot(412)
plot(WKLT2(2,:)) ;
title('2nd basis (eigen)vector')

subplot(413)
plot(WKLT2(3,:)) ;
title('3rd basis (eigen)vector')

subplot(414)
plot(WKLT2(5,:)) ;
title('4th basis (eigen)vector')

% %%%%%%%%%%%%%%%%%%%%%%%%

figure(5)
subplot(311)
plot(s3)
title('SIGNAL s3 - whole')

subplot(312)
plot(SS3(7,:)) ;
title('SIGNAL s3 - whole')

subplot(313)
stem(diag(D3)) ;
title('SIGNAL s3 - eigenvalues of cov.matrix')

figure(6)
subplot(411)
plot(WKLT3(1,:))
title('SIGNAL s3 - 1st basis (eigen)vector (most significant)')

subplot(412)
plot(WKLT3(2,:)) ;
title('2nd basis (eigen)vector')

subplot(413)
plot(WKLT3(3,:)) ;
title('3rd basis (eigen)vector')

subplot(414)
plot(WKLT3(4,:)) ;
title('4th basis (eigen)vector')


% Frame for transformation & transforms

start=1070;
frame1=s1(start+1:start+wlen);
frame2=s2(start+1:start+wlen);
frame3=s3(start+1:start+wlen);

Y_KLT1= WKLT1*frame1 ;
Y_KLT2= WKLT2*frame2 ;
Y_KLT3= WKLT3*frame3 ;

Y_KLT21= WKLT1*frame2 ;
Y_KLT23= WKLT3*frame2 ;

Y_DCT1 = dct(frame1) ;
Y_DCT2 = dct(frame2) ;
Y_DCT3 = dct(frame3) ;

Y_FFT1 = fft(frame1) ;
Y_FFT2 = fft(frame2) ;
Y_FFT3 = fft(frame3) ;


figure(7)
subplot(221)
plot(frame1)
title('frame1 (from signal s1)')

subplot(223)
plot(Y_KLT1,'.-') ;
title('KLT transform of frame1')

subplot(222)
plot(abs(Y_FFT1),'.-')
title('Magnitude FFT of frame1')

subplot(224)
plot(Y_DCT1,'.-') ;
title('DCT transform of frame1')

figure(8)
subplot(221)
plot(frame2)
title('frame2 (from signal s2)')

subplot(223)
plot(Y_KLT2,'.-') ;
title('KLT transform of frame2')

subplot(222)
plot(abs(Y_FFT2),'.-')
title('Magnitude FFT of frame2')

subplot(224)
plot(Y_DCT2,'.-') ;
title('DCT transform of frame2')

figure(9)
subplot(221)
plot(frame3)
title('frame3 (from signal s3)')

subplot(223)
plot(Y_KLT3,'.-') ;
title('KLT transform of frame3')

subplot(222)
plot(abs(Y_FFT3),'.-')
title('Magnitude FFT of frame3')

subplot(224)
plot(Y_DCT3,'.-') ;
title('DCT transform of frame3')

figure(10)
subplot(311)
plot(Y_KLT21,'.-')
title('KLT transform of frame2 with basis WKLT1 (BAD)')

subplot(312)
plot(Y_KLT2,'.-') ;
title('KLT transform of frame2 with basis WKLT2 (OK)')

subplot(313)
plot(Y_KLT23,'.-')
title('KLT transform of frame2 with basis WKLT3 (BAD)')

% INVERSE TRANSFORM

IWKLT1 = inv(WKLT1) ;
IWKLT2 = inv(WKLT2) ;
IWKLT3 = inv(WKLT3) ;

outframe1=IWKLT1*Y_KLT1 ;
outframecomps1=[];
for i=1:6,
    mask=[ones(i,1); zeros(wlen-i,1)];
    outframecomps1=[outframecomps1 IWKLT1*(Y_KLT1.*mask)] ;
end

outframe2=IWKLT2*Y_KLT2 ;
outframecomps2=[];
for i=1:6,
    mask=[ones(i,1); zeros(wlen-i,1)];
    outframecomps2=[outframecomps2 IWKLT2*(Y_KLT2.*mask)] ;
end

outframe3=IWKLT3*Y_KLT3 ;
outframecomps3=[];
for i=1:6,
    mask=[ones(i,1); zeros(wlen-i,1)];
    outframecomps3=[outframecomps3 IWKLT3*(Y_KLT3.*mask)] ;
end

figure(11)
subplot(211)
plot(frame1)
hold on
plot(outframe1,'r--')
plot(outframecomps1(:,6),'k-.')
hold off
title('frame1 - original & inverse from all (red) / 6 (black)')
subplot(212)
plot(outframecomps1)
title('inverse transform from 1-6 components')


figure(12)
subplot(211)
plot(frame2)
hold on
plot(outframe2,'r--')
plot(outframecomps2(:,6),'k-.')
hold off
title('frame2 - original & inverse from all (red) / 6 (black)')
subplot(212)
plot(outframecomps2)
title('inverse transform from 1-6 components')

figure(13)
subplot(211)
plot(frame3)
hold on
plot(outframe3,'r--')
plot(outframecomps3(:,6),'k-.')
hold off
title('frame3 - original & inverse from all (red) / 6 (black)')
subplot(212)
plot(outframecomps3)
title('inverse transform from 1-6 components')

%% COMPRESSION USING FFT, DCT and KLT

% frame1 from sig1
% DFT compression
i=0;
Pf=mean(frame1.^2);
Pout=0;
while ( Pout/Pf < 0.95 ),
    i=i+1;
    mask=[ones(i,1); zeros(wlen-2*i+1,1); ones(i-1,1);];
    outframefft1=ifft(Y_FFT1.*mask) ;
    Pout=mean(outframefft1.^2);
end

figure(21)
subplot(311)
plot(frame1)
hold on
plot(outframefft1,'k-.')
hold off
title(['FFT compression: components=' num2str(2*i-1) ', ' num2str(100*Pout/Pf) '% of power'])

% DCT compression
i=0;
Pf=mean(frame1.^2);
Pout=0;
while ( Pout/Pf < 0.95 ),
    i=i+1;
    mask=[ones(i,1); zeros(wlen-i,1)];
    outframedct1=idct(Y_DCT1.*mask) ;
    Pout=mean(outframedct1.^2);
end

subplot(312)
plot(frame1)
hold on
plot(outframedct1,'k-.')
hold off
title(['DCT compression: components=' num2str(i) ', ' num2str(100*Pout/Pf) '% of power'])

% KLT compression
i=0;
Pf=mean(frame1.^2);
Pout=0;
while ( Pout/Pf < 0.95 ),
    i=i+1;
    mask=[ones(i,1); zeros(wlen-i,1)];
    outframeklt1=IWKLT1*(Y_KLT1.*mask) ;
    Pout=mean(outframeklt1.^2);
end

subplot(313)
plot(frame1)
hold on
plot(outframeklt1,'k-.')
hold off
title(['KLT compression: components=' num2str(i) ', ' num2str(100*Pout/Pf) '% of power'])


% frame2 from sig2
% DFT compression
i=0;
Pf=mean(frame2.^2);
Pout=0;
while ( Pout/Pf < 0.95 ),
    i=i+1;
    mask=[ones(i,1); zeros(wlen-2*i+1,1); ones(i-1,1);];
    outframefft2=ifft(Y_FFT2.*mask) ;
    Pout=mean(outframefft2.^2);
end

figure(22)
subplot(311)
plot(frame2)
hold on
plot(outframefft2,'k-.')
hold off
title(['FFT compression: components=' num2str(2*i-1) ', ' num2str(100*Pout/Pf) '% of power'])

% DCT compression
i=0;
Pf=mean(frame2.^2);
Pout=0;
while ( Pout/Pf < 0.95 ),
    i=i+1;
    mask=[ones(i,1); zeros(wlen-i,1)];
    outframedct2=idct(Y_DCT2.*mask) ;
    Pout=mean(outframedct2.^2);
end

subplot(312)
plot(frame2)
hold on
plot(outframedct2,'k-.')
hold off
title(['DCT compression: components=' num2str(i) ', ' num2str(100*Pout/Pf) '% of power'])

% KLT compression
i=0;
Pf=mean(frame2.^2);
Pout=0;
while ( Pout/Pf < 0.95 ),
    i=i+1;
    mask=[ones(i,1); zeros(wlen-i,1)];
    outframeklt2=IWKLT2*(Y_KLT2.*mask) ;
    Pout=mean(outframeklt2.^2);
end

subplot(313)
plot(frame2)
hold on
plot(outframeklt2,'k-.')
hold off
title(['KLT compression: components=' num2str(i) ', ' num2str(100*Pout/Pf) '% of power'])


% frame3 from sig3
% DFT compression
i=0;
Pf=mean(frame3.^2);
Pout=0;
while ( Pout/Pf < 0.95 ),
    i=i+1;
    mask=[ones(i,1); zeros(wlen-2*i+1,1); ones(i-1,1);];
    outframefft3=ifft(Y_FFT3.*mask) ;
    Pout=mean(outframefft3.^2);
end

figure(23)
subplot(311)
plot(frame3)
hold on
plot(outframefft3,'k-.')
hold off
title(['FFT compression: components=' num2str(2*i-1) ', ' num2str(100*Pout/Pf) '% of power'])

% DCT compression
i=0;
Pf=mean(frame3.^2);
Pout=0;
while ( Pout/Pf < 0.95 ),
    i=i+1;
    mask=[ones(i,1); zeros(wlen-i,1)];
    outframedct3=idct(Y_DCT3.*mask) ;
    Pout=mean(outframedct3.^2);
end

subplot(312)
plot(frame3)
hold on
plot(outframedct3,'k-.')
hold off
title(['DCT compression: components=' num2str(i) ', ' num2str(100*Pout/Pf) '% of power'])

% KLT compression
i=0;
Pf=mean(frame3.^2);
Pout=0;
while ( Pout/Pf < 0.95 ),
    i=i+1;
    mask=[ones(i,1); zeros(wlen-i,1)];
    outframeklt3=IWKLT3*(Y_KLT3.*mask) ;
    Pout=mean(outframeklt3.^2);
end

subplot(313)
plot(frame3)
hold on
plot(outframeklt3,'k-.')
hold off
title(['KLT compression: components=' num2str(i) ', ' num2str(100*Pout/Pf) '% of power'])



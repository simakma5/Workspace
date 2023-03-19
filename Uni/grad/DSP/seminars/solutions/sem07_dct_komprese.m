% Cviceni na Diskretni Kosinovou Transformaci
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KOMPRESE POMOCI DFT a DCT

% clear all;
% close all ;

x=loadbin('frame.bin');
N=length(x);

N=256;
n=0:N-1;
x=0.99.^n.*cos(0.02*pi*n);

x=x(:);

figure(1)
subplot(211)
plot(x)
aa=axis ; axis([0 N aa(3) aa(4)]) ;


XC1=dctxc1(x);
%XC1=sqrt(N/2)*dct(x);

subplot(212)
semilogy(1:N,abs(XC1),'r') ;
aa=axis ; axis([0 N 1e-2 1e2]) ;

X=fft(x); 
hold on
semilogy(1:N/2,abs(X(1:N/2)),'b'); 
hold off
aa=axis ; axis([0 N/2 1e-2 1e2]) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATION OF ENERGY IN STARTING BLOCK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)

pwrdct=[];
pwrdft=[];

for k=1:N  ;
  pwrdct=[pwrdct ; sum(XC1(1:k).^2)];
end

for k=1:N/2  ;
  pwrdft=[pwrdft ; sum(abs(X(1:k)).^2)];
end

plot(1:N,pwrdct,'r');
hold on
plot(1:2:2*k,pwrdft,'b');
hold off
xlabel('The first k-componetns');
ylabel('Energy in block');
title('red - DCT based power compression, blue - DFT based power compression')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRUNCATION OF END LOW ENERGY BLOCK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)

errdct=[];
errdft=[];

for k=0:N/2 -1 ;

  XC=[XC1(1:end-k); zeros(k,1)];
  xx=idctxc1(XC) ;
  %% xx=idct(XC);
  
  errdct=[errdct ; mean((x-xx).^2)];
  
  XX=X ; 
  XX(N/2+1)=0;
  XX(N/2+2:N/2+1+k)=zeros(k,1);
  XX(N/2-k+1:N/2)=zeros(k,1);
  
  xxx=ifft(XX);
  
  errdft=[errdft ; mean((x-xxx).^2)];
  
end

semilogy(0:k,errdct,'r',0:k,errdft,'b');
xlabel('Cutting of k-componetns');
ylabel('mean square error');
title('red - DCT based compression, blue - DFT based compression')


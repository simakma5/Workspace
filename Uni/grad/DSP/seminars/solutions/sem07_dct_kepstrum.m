% Cviceni na Diskretni Kosinovou Transformaci
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % PODCAST - Vypocet realneho kepstra pomoci DCT
% 
frame=loadbin('frame.bin');

c1=rceps(frame);
N=length(frame);

figure(1)

subplot(211)
plot(c1(1:N/2))
title('Realne kepstrum pomoci rceps')

X=log(abs(fft(frame)));
c2=real(ifft(X));

k=1:N/2+1;
X=X(k);
c3=idctxc1(X) ;



subplot(212);
plot(c2(1:N/2))
hold on
plot([c3],'r')
hold off
title('Realne kepstrum pomoci idft (blue) a idct (red) - podle idctx1.m (2N-2 extense)')



figure(2)

subplot(211)
plot(c1(1:N/2))
title('Realne kepstrum pomoci rceps')

X=log(abs(fft(frame)));
c2=real(ifft(X));

k=1:N/2+1;
X=X(k);
c3=1/sqrt(N-1)*idct(X);

subplot(212);
plot(c2(1:N/2))
hold on
plot([c3],'r')
hold off
title('Realne kepstrum pomoci idft (blue) a idct (red) - idct v MATLABu (2N extense)')


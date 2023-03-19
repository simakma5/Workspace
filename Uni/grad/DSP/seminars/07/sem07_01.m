%%
clear, clc, close all;

%% Task 1 using DCT-1
x = [10 7 4 1 -2 -5 -8 -5 -2 1];
N = length(x);
n = 0:N-1;

x1 = [x x(end-1:-1:2)];
% x1 = [x fliplr(x(2:end))];  % should be the same thing
n1 = 0:2*N-2-1;

Xdct = dctxc1(x);
Xfft = real(fft(x1));

figure(1)
subplot(221);
stem(n,x)
xlim([0 2*N-1])
title('Original sequence x')

subplot(222);
plot(n,Xdct)
xlim([0 2*N-1])
title('dctxc1(x)')

subplot(223);
stem(n1,x1)
xlim([0 2*N-1])
title('Sequence x1')

subplot(224);
plot(n1,Xfft)
xlim([0 2*N-1])
title('fft(x1)')

%% Task 1 using DCT-2
x = [10 7 4 1 -2 -5 -8 -5 -2 1];
N = length(x);
n = 0:N-1;

x1 = [x x(end:-1:1)];
% x1 = [x fliplr(x(1:end))];  % should be the same thing
n1 = 0:2*N-1;

Xdct = dctxc2(x);
Xfft = real(fft(x1).*exp(1i*pi*n1/(2*N)));

figure(2)
subplot(221);
stem(n,x)
xlim([0 2*N-1])
title('Original sequence x')

subplot(222);
plot(n,Xdct)
xlim([0 2*N-1])
title('dctxc2(x)')

subplot(223);
stem(n1,x1)
xlim([0 2*N-1])
title('Sequence x1')

subplot(224);
plot(n1,Xfft)
xlim([0 2*N-1])
title('fft(x1)')

figure(3)
plot(n,Xdct)
xlim([0 2*N-1])
hold on
plot(n1,Xfft)
xlim([0 2*N-1])
hold off
% nesedi :(


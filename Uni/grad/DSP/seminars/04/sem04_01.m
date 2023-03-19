% we'll want to create a channel communication, where we have signal
% which is then split into two channels - the offbranched signal is then
% delayed by a signal block
fs = 16e03;
s = loadbin('SA001S01.CS0');

delay = 20;
sim1s1 = s(delay+1:end);
sim1s2 = s(1:end-delay);

figure(1)
subplot(211);
plot([sim1s1 sim1s2])

G = 0.001;
M = 30;  % order of the filter = 1/2 of the delay
Wn = 0.05;

b = fir1(M,Wn,'high');

sim2s1 = s + G*randn(length(s),1);
sim2s2 = filter(b,1,s) + G*randn(length(s),1);

subplot(212);
plot([sim2s1 sim2s2])

% inspect the frequency response of the created FIR filter
figure(2)
freqz(b,1)

% check impulse response (matches the paper notes)
figure(3)
impz(b,1)

%%

% choose a segment of the signal (the whole is uninteresting)
siglen = 3e3;
start = 2e4;

x1 = sim1s1(start+1:start+siglen);
x2 = sim1s2(start+1:start+siglen);
% x1 = sim2s1(start+1:start+siglen);
% x2 = sim2s2(start+1:start+siglen);

[R12,lags] = xcorr(x1,x2);

[~,ix] = max(R12);  % suppress first argument just to retreve the point
delay = lags(ix);

figure(4)
plot(lags,R12), grid
title(['delay = ' num2str(delay) ' samples'])
xline(delay, 'r')

%%

wlen = 512;
noverlap = wlen/2;
N = wlen;

% cross-PSD determines delay and shit
[S12, ff] = cpsd(x1,x2,wlen,noverlap,N,fs);

figure(5)
subplot(211);
plot(ff,10*log10(abs(S12)))

subplot(212);
plot(ff,unwrap(angle(S12)))

% finish at home: calculate delay using the slope of the angle plot
%  - might need to set proper scaling: slope*(2pi/fs)??

%%
wlen2 = 128;
noverlap2 = 3*wlen2/4;
nfft2 = max([wlen2,256]);
[S12b, ffb] = cpsd(x1,x2,wlen2,noverlap2,nfft2,fs);

fmax = 3500;
nfftmax = find(ffb>=fmax,1);



d2 = mean(phd(1:nfftmax));
title(['delay = ' num2str(d2) 'samples'])

figure(11)
subplot(211);
plot([x1 x2])

subplot(212);
mscogere(x1,x2,wlen2,noverlap2,nfft2,fs)
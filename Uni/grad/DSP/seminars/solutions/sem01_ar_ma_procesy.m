%   Demo AR / MA processes (signal modeling)

% Ad 1) typ DP - proces 1. radu

% fs=8000;
% Bo=300 ;
% 
% zo=exp(-Bo*pi/fs);
% a=[1 -zo];

% 1st order models

fs=8000;
Bo=300 ;
zo=exp(-Bo*pi/fs);

% 1st order low-frequency AR model
b=1 ;
a=[1 -zo];
type='1st order AR - LF';

% 1st order high-frequency AR model
% b=1 ;
% a=[1 zo];
% type='1st order AR - HF';

% 1st order low-frequency MA model
% a=1 ;
% b=[1 zo];
% type='1st order MA - LF';

% 1st order high-frequency MA model
% a=1 ;
% b=[1 -zo];
% type='1st order MA - HF';

% 2nd order model AR model - band-pass (rezonator)
% p1 = 0.9 * exp ( j * 0.9 ) ;
% p2 = 0.9 * exp ( - j * 0.9 ) ;
% a=poly([p1 p2]);
% b=1 ;
% type='2nd order AR - BP';


figure(1)
freqz(b,a,1000,fs)
subplot(211)
title(['Frequency reponse - ' type])

figure(2)
zplane(b,a)

% signal modeling (filtering of white noise)

N=10000;
no=randn(N,1);

nc=filter(b,a,no);

wlen=512;

figure(3)
spectrogram(nc,wlen,[],[],fs,'yaxis')
title('DFT spectrogram');
colorbar off;
colormap jet;

figure(4)
NC=fft(nc(1:wlen).*hamming(wlen),wlen);
ff=(0:wlen-1)./wlen*fs;
plot(ff,10*log10((abs(NC).^2)./wlen));
title('short-time power spectrum in dBs (periodogram)');

figure(5)
pwelch(nc,wlen,[],[],fs);
title('PSD in dB using Welch method');




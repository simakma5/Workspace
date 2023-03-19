% Multi-band signal processing - filter banks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design of filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs=16000;

% MATLAB prototype
%%%%%%%%%%%%%%%%%%%%%%%%
M=50;
b0=fir1(M,0.5,'low');

% Perfect reconstruction filter bank (QMF filters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% b0= [ -0.0059, 0.0129, 0.0013, -0.0274, 0.0086, 0.0510, -0.0338, -0.1001, 0.1243, 0.4688, 0.4688, 0.1243, -0.1001, -0.0338, 0.0510, 0.0086, -0.0274, 0.0013, 0.0129, -0.0059 ] ;
% M=length(b0)-1;
    
H0=freqz(b0,1,1000) ;
H0abs=20*log10(abs(H0));

b1 = b0.*((-1).^(0:M));
H1=freqz(b1,1,1000) ;
H1abs=20*log10(abs(H1));

bg0 = fliplr(b0);
ag0 = [zeros(1,M) 1] ;
G0=freqz(bg0,ag0,1000);
G0abs=20*log10(abs(G0));

bg1 = fliplr(b1);
ag1 = [zeros(1,M) 1] ;
G1=freqz(bg1,ag1,1000);
G1abs=20*log10(abs(G1));

figure(1)
plot(H0abs,'b');
grid
hold on
plot(H1abs,'r');
plot(G0abs,'g--');
plot(G1abs,'k--');
hold off
title('Magnitude responses of filters H0, H1, G0, G1')

figure(2)
plot(unwrap(angle(H0)),'b');
grid
hold on
plot(unwrap(angle(H1)),'r');
plot(unwrap(angle(G0)),'g--');
plot(unwrap(angle(G1)),'k--');
hold off
title('Phase responses of filters H0, H1, G0, G1')


figure(3)
subplot(221)
stem(b0);
title('Impulse response of H0(z)')
subplot(223)
stem(b1)
title('Impulse response of H1(z)')
subplot(222)
stem(bg0);
title('Impulse response of G0(z)')
subplot(224)
stem(bg1)
title('Impulse response of G1(z)')

figure(4)
subplot(211)
plot(abs(H0).^2+abs(H1).^2);
title('|H0|^2 + |H1|^2')
a=axis;
subplot(212)
plot(abs(H0)+abs(H1));
title('|H0| + |H1|')
axis(a)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtrace do pasem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s=loadbin('s0001.bin');

s0=filter(b0,1,s);
s1=filter(b1,1,s);

wlen=floor(0.032*fs);

figure(10)
subplot(211);
tt=(1:length(s))./fs;
plot(tt,s); xlim([0 max(tt)])
title('Input signal and its spectrogram');
subplot(212);
spectrogram(s,wlen,[],[],fs,'yaxis'); 
colorbar off;



figure(11)
subplot(211);
spectrogram(s0,wlen,[],[],fs,'yaxis'); 
colorbar off;
title('Signal s0 (after filtering to LFB)');
subplot(212);
spectrogram(s1,wlen,[],[],fs,'yaxis'); 
colorbar off;
title('Signal s1 (after filtering to HFB)');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z0=s0(1:2:length(s0)) ;
z1=s1(1:2:length(s1)) ;

figure(12)
subplot(211);
spectrogram(z0,wlen/2,[],[],fs/2,'yaxis'); 
colorbar off;
title('Spectrogram of decimated signal in LFB (z0)');
subplot(212);
spectrogram(z1,wlen/2,[],[],fs/2,'yaxis'); 
colorbar off;
title('Spectrogram of decimated signal in HFB (z1)');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolation - zero interlacing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0=zeros(2*length(z0),1) ;
x1=zeros(2*length(z1),1) ;

x0(1:2:end)=z0 ;
x1(1:2:end)=z1 ;

figure(13)
subplot(211);
spectrogram(x0,wlen,[],[],fs,'yaxis'); 
colorbar off;
title('Spectrogram of zero-interpolated signal in LFB (x0)');
subplot(212);
spectrogram(x1,wlen,[],[],fs,'yaxis'); 
colorbar off;
title('Spectrogram of zero-interpolated signal in HFB (x1)');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. krok interpolace - filtrace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% y0=filter(bg0,ag0,x0) ;
% y1=filter(bg1,ag1,x1) ;
y0=filter(bg0,1,x0) ;
y1=filter(bg1,1,x1) ;

figure(14)
subplot(211);
spectrogram(y0,wlen,[],[],fs,'yaxis'); 
colorbar off;
title('Spectrogram filtered interpolated signal in LFB');
subplot(212);
spectrogram(y1,wlen,[],[],fs,'yaxis'); 
colorbar off;
title('Spectrogram filtered interpolated signal in HFB');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Secteni signalu z obou pasem v puvodni fs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y=2*(y0+y1) ;
y=y(1:length(s));

figure(15)
subplot(211);
tt=(1:length(y))./fs ;
plot(tt,y); xlim([0 max(tt)]) ;
title('Output signal');
subplot(212);
spectrogram(y,wlen,[],[],fs,'yaxis'); 
colorbar off;
title('Spectrogram of output signal');

%% ERROR SIGNAL ANALYSIS

e=s(1:end-M)-y(M+1:end);
% e=s-y;
figure(20)
% plot(e(11000:13000))
plot(e)

Pe=mean(e.^2);
Ps=mean(s.^2);
SNR_e=10*log10(Ps/Pe);

title(['Error signal, SNR_e = ' num2str(SNR_e) ' dB'])





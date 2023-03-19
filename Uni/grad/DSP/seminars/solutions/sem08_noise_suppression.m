% Various methods of noise suppresion

% Loading of harmonic components
load sink.mat
fs=16000;
ttt=(1:4000)./fs;
% harmonic components
clean=sink(1,:)+sink(2,:)+sink(3,:);
% % %% More harmonic components
% clean=sum(sink) ;

% clean=loadbin('vm0.bin');  fs=16000 ;


clean=clean(:);
slen=length(clean);

% noi=loadbin('nc1.bin');
noi=loadbin('nc2.bin');
% noi=loadbin('nc3.bin');

% Addition of noise
% %%%%%%%%%%%%%%%%%

scale=0.3;
noi=noi(:);
noi=scale*noi(1:slen);

sig=clean+noi;


% Loading of noisy speech
% sig=loadbin('SA110992_auto1.CS0');  fs=16000 ;
% sig=loadbin('SA107S06_auto2.CS0');  fs=16000 ;
% sig=loadbin('ma034014_auto3.ils');  fs=16000 ;


% Vypocet poctu oken
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
slen=length(sig);
wlen=512 ;
wstep=wlen/2 ;
wnum=(slen-wlen)/wstep+1 ;

w=hamming(wlen);
% w=rectwin(wlen);
% w=hanning(wlen);
% w=blackman(wlen);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation of magnitude spectrum of noise (for SS)
% Estimation of PSD of noise (for WF)
%   z prvniho casti signalu bez reci
%   (prvnich 30 segmentu, nebo wnum pro kratsi signaly)
winit =min([wnum 30]) ;
Nabs=[] ;
Nabs2=[] ;
Sabs2=[] ;
Ndct=[];
Sdct=[];
for i=1:winit,
  ii=(i-1)*wstep+1;
  jj=(i-1)*wstep+wlen;
  nframe=noi(ii:jj).*w;
  N=fft(nframe);
  Nabs = [Nabs abs(N)] ;
  Nabs2 = [Nabs2 (abs(N).^2)./wlen] ;
  Ndct = [Ndct dct(nframe)];
  
  sframe=clean(ii:jj).*w;
  S=fft(sframe);
  Sabs2 = [Sabs2 (abs(S).^2)./wlen] ;
  Sdct = [Sdct dct(sframe)];
end;

% For Spectral Subtraction
Navg=mean(Nabs,2);

% For Wiener filtering
% Npsd=mean(Nabs2,2);    
% Spsd=mean(Sabs2,2);
     % Estimation with direct averaging
Npsd=pwelch(noi(1:wstep*winit),wlen,wstep,wlen,'twosided')/mean(hamming(wlen).^2);
Spsd=pwelch(clean(1:wstep*winit),wlen,wstep,wlen,'twosided')/mean(hamming(wlen).^2);
     % Estimation using pwelch (ATTENTION! Normalization to power of used
     % weighting window must be applied
H=sqrt(Spsd./(Spsd+Npsd));
H=(Spsd./(Spsd+Npsd));

% For DCT based noise suppression
Ndct2avg=mean(Ndct.^2,2);
Sdct2avg=mean(Sdct.^2,2);
aprio_snr=(Sdct2avg)./(Ndct2avg);
G=sqrt(aprio_snr./(aprio_snr+1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generovani nulove vystupni posloupnosti
out = zeros(slen,1);
out2 = zeros(slen,1);
out3 = zeros(slen,1);

% Main cycle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:wnum,

  % short-time frame selection
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ii=(i-1)*wstep+1;
  jj=(i-1)*wstep+wlen;

  frame=sig(ii:jj).*w;
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % FFT-based spectral subtraction
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  X=fft(frame);
  Xabs = abs(X) ;
  Xangle = angle(X) ;
  
  Yabs = Xabs - Navg ;
  Yabs = ( Yabs + abs(Yabs) )./ 2 ;   % Half-wave rectification
  %  Yabs = abs( Yabs ) ;             % Full-wave rectification
  
  yframe = real ( ifft( Yabs .* exp (j*Xangle) ) ) ;
  % Addition to output signal (implementation of OLA with general window)
  out(ii:jj)=out(ii:jj)+yframe ;
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Wiener filter (basic)
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  y2frame = real ( ifft( X .* H ) ) ;
  out2(ii:jj)=out2(ii:jj)+y2frame ;
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % DCT-based noise suppression
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  X=dct(frame);

  Y=G.*X ;
  
  y3frame = idct ( Y ) ;
  % Addition to output signal (implementation of OLA with general window)
  out3(ii:jj)=out3(ii:jj)+y3frame ;

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalization of output signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scale2 = sum ( w(1:wstep:wlen) ) ;
out = out ./ scale2 ;
out2 = out2 ./ scale2 ;
out3 = out3 ./ scale2 ;

Pclean = mean(clean(wlen:end-wlen).^2);
Pnoi = mean(noi(wlen:end-wlen).^2);

SNRin = 10*log10(Pclean/Pnoi);

noiout=out-clean;
noiout2=out2-clean;
noiout3=out3-clean;

Pnoiout=mean(noiout(wlen:end-wlen).^2);
Pnoiout2=mean(noiout2(wlen:end-wlen).^2);
Pnoiout3=mean(noiout3(wlen:end-wlen).^2);

SNRout=10*log10(Pclean/Pnoiout);
SNRout2=10*log10(Pclean/Pnoiout2);
SNRout3=10*log10(Pclean/Pnoiout3);


% Graphs for Spectral subtraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf ;
subplot(211);
plot(sig);
title(['Input signal: SNRin = ' num2str(SNRin) ' dB']);

subplot(212);
plot(out) ;
title(['Output SSub:  SNRout = ' num2str(SNRout) ' dB; SNRE = ' num2str(SNRout-SNRin) ' dB']);

figure(2)
subplot(211);
spectrogram(sig,wlen,[],[],fs,'yaxis'); 
title('Input signal');
colorbar off

subplot(212);
spectrogram(out,wlen,[],[],fs,'yaxis'); 
title(['Output SSub: SNRE = ' num2str(SNRout-SNRin) ' dB']);
colorbar off

% Graphs for Wiener filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3); clf ;
subplot(211);
plot(sig);
title(['Input signal: SNRin = ' num2str(SNRin) ' dB']);

subplot(212);
plot(out2) ;
title(['Output WF: SNRout2 = ' num2str(SNRout2) ' dB; SNRE = ' num2str(SNRout2-SNRin) ' dB']);

figure(4)
subplot(211);
spectrogram(sig,wlen,[],[],fs,'yaxis'); 
title('Input signal');
colorbar off

subplot(212);
spectrogram(out2,wlen,[],[],fs,'yaxis'); 
title(['Output WF: SNRE = ' num2str(SNRout2-SNRin) ' dB']);
colorbar off

% Graphs for Spectral subtraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5); clf ;
subplot(211);
plot(sig);
title(['Input signal: SNRin = ' num2str(SNRin) ' dB']);

subplot(212);
plot(out3) ;
title(['Output DCT: SNRout3 = ' num2str(SNRout3) ' dB; SNRE = ' num2str(SNRout3-SNRin) ' dB']);

figure(6)
subplot(211);
spectrogram(sig,wlen,[],[],fs,'yaxis'); 
title('Input signal');
colorbar off

subplot(212);
spectrogram(out3,wlen,[],[],fs,'yaxis'); 
title(['Output DCT: SNRE = ' num2str(SNRout3-SNRin) ' dB']);
colorbar off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse filtering

sig=loadbin('SA001S01.CS0');   speechpart=1.8e4:8e4 ;
fs=16000;
slen=length(sig);

% Vypocet poctu oken
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
slen=length(sig);
wlen=512 ;
wstep=wlen/2 ;
wnum=(slen-wlen)/wstep+1 ;

snrest=0;
% snrest=1;

if snrest==1,

    % SNR estimation Cycle for power comuptation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:wnum,
        % short-time frame selection
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ii=(i-1)*wstep+1;
        jj=(i-1)*wstep+wlen;
        frame=sig(ii:jj);
        Pvdb(i)=10*log10(mean(frame.^2)) ;
    end
    Pvdbsort=sort(Pvdb);  estpart=round(0.1*length(Pvdbsort));
    SNRorig = mean(Pvdbsort(end-estpart:end))-mean(Pvdbsort(1:estpart));
 
else
    SNRorig=-inf;
end

% %%% Preemphasis
% hb=[1 -0.95];
% ha=1;

% %%% Band channel distortion
pabs=0.98;
pphi=11*pi/12;
p1 = pabs * exp ( j * pphi ) ;
p2 = pabs * exp ( - j * pphi ) ;
ha=poly([p1 p2]);
hb=1;



figure(1)
freqz(hb,ha,1000)
title('Frequency response of CHANNEL MODEL')

figure(2)
freqz(ha,hb,1000)
title('Frequency response of INVERSE FILTER')

sigout=filter(hb,ha,sig);   % Channel filtering
Pout=mean(sigout(speechpart).^2);

scale=0.002;
% scale=0;
chnoise=scale*randn(slen,1);
sigoutnoise=sigout+chnoise;

Pnoi=mean(chnoise.^2);
SNR=10*log10(Pout/Pnoi);

siginv=filter(ha,hb,sigoutnoise);  % Inverse filtering

if snrest==1,
    
    % SNR estimation Cycle for power comuptation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:wnum,
        % short-time frame selection
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ii=(i-1)*wstep+1;
        jj=(i-1)*wstep+wlen;
        frame=siginv(ii:jj);
        Pvdb(i)=10*log10(mean(frame.^2)) ;
    end
    Pvdbsort=sort(Pvdb);  estpart=round(0.1*length(Pvdbsort));
    SNRinv = mean(Pvdbsort(end-estpart:end))-mean(Pvdbsort(1:estpart));
 
else

    Psig=mean(sig(speechpart).^2);
    Pinvnoi=mean((siginv-sig).^2);
    SNRinv=10*log10(Psig/Pinvnoi);
    
end

figure(11)
subplot(211)
plot(sig) ;
title(['Input signal']);

subplot(212)
spectrogram(sig,512,[],[],fs,'yaxis');
colorbar off
title(['SNR_{orig}=' num2str(SNRorig) ' dB']);

figure(12)
subplot(211)
plot(sigout) ;
title(['Channel signal']);

subplot(212)
spectrogram(sigout,512,[],[],fs,'yaxis');
colorbar off
title('Spectrogram of channel signal')

figure(13)
subplot(211)
plot(sigoutnoise) ;
title(['Channel signal with noise']);

subplot(212)
spectrogram(sigoutnoise,512,[],[],fs,'yaxis');
colorbar off
title(['SNR = ' num2str(SNR) ' dB'])


figure(14)
subplot(211)
plot(siginv) ;
title(['Inverse filtered signal']);

subplot(212)
spectrogram(siginv,512,[],[],fs,'yaxis');
colorbar off
title(['SNR_{inv} = ' num2str(SNRinv) ' dB'])


% Listening
% sound([sig; sigout; sigoutnoise; siginv],fs)

%% Wienerova fitrace
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computation of the number of short-time frames 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
slen=length(sig);
wlen=512 ;
wstep=wlen/2 ;
wnum=(slen-wlen)/wstep+1 ;

w=hamming(wlen);
% w=rectwin(wlen);
% w=hanning(wlen);
% w=blackman(wlen);


H=freqz(hb,ha,wlen,'whole');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation of PSD of noise (for WF)
%   from the first part of the signal without speech
%   (the first 30 frames)


winit =min([wnum 30]) ;
Nabs=[] ;
Nabs2=[] ;
Sabs2=[] ;
for i=1:winit,
   ii=(i-1)*wstep+1;
   jj=(i-1)*wstep+wlen;
   nframe=sigoutnoise(ii:jj).*w;
   N=fft(nframe);
   Nabs2 = [Nabs2 (abs(N).^2)./wlen] ;
end;
% For Wiener filtering
Npsd=mean(Nabs2,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creation of initial output sequence with zeros
outfreq = zeros(slen,1);
outwf = zeros(slen,1);

% Main cycle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:wnum,

  % short-time frame selection
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ii=(i-1)*wstep+1;
  jj=(i-1)*wstep+wlen;

  frame=sigoutnoise(ii:jj).*w;
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % computation of short-time frame spectrum
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  X=fft(frame);
  Xpsd=(abs(X).^2)./wlen;
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % pure inverse filtering in frequency domain
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Y = 1./H.*X ;
  yframe = real ( ifft( Y ) ) ;
  outfreq(ii:jj)=outfreq(ii:jj)+yframe ;
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % inverse filtering with Wiener filter (with distortion model)
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Inverse filtering with WF
    %     Npsd .... input noise PSD estimation
    %     Xpsd .... power spectrum of current short-time frame

  SNRprio = (Xpsd-Npsd)./Npsd ;
  SNRprio(SNRprio<0) = 0 ;
  Hk = 1./H.*(abs(H).^2)./(abs(H).^2 + 1./SNRprio) ;
    
  Y = Hk.*X ;
  
  yframe = real ( ifft( Y ) ) ;
  outwf(ii:jj)=outwf(ii:jj)+yframe ;
  
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalization of output signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scale2 = sum ( w(1:wstep:wlen) ) ;
outfreq = outfreq ./ scale2 ;
outwf = outwf ./ scale2 ;

if snrest==1,
    
    % SNR estimation Cycle for power comuptation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:wnum,
        % short-time frame selection
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ii=(i-1)*wstep+1;
        jj=(i-1)*wstep+wlen;
        frame1=outfreq(ii:jj);
        frame2=outwf(ii:jj);
        Pvdbfreq(i)=10*log10(mean(frame1.^2)) ;
        Pvdbwf(i)=10*log10(mean(frame2.^2)) ;
    end
    Pvdbsort=sort(Pvdbfreq);  estpart=round(0.1*length(Pvdbsort));
    SNRinvfreq = mean(Pvdbsort(end-estpart:end))-mean(Pvdbsort(1:estpart));
    Pvdbwf=sort(Pvdbwf);  estpart=round(0.1*length(Pvdbsort));
    SNRinvwf = mean(Pvdbsort(end-estpart:end))-mean(Pvdbsort(1:estpart));
 
else

    Pninvfreq=mean((outfreq-sig).^2);
    SNRinvfreq=10*log10(Psig/Pninvfreq);

    Pninvwf=mean((outwf-sig).^2);
    SNRinvwf=10*log10(Psig/Pninvwf);
    
end


figure(15)
subplot(211)
plot(outfreq) ;
title(['Output signal - inverse in frequency domain']);

subplot(212)
spectrogram(outfreq,512,[],[],fs,'yaxis');
colorbar off
title(['SNR_{inv,freq} = ' num2str(SNRinvfreq) ' dB'])

figure(16)
subplot(211)
plot(outwf) ;
title(['Output signal - WF (with signal model)']);

subplot(212)
spectrogram(outwf,512,[],[],fs,'yaxis');
colorbar off
title(['SNR_{inv,wf}  = ' num2str(SNRinvwf) ' dB'])



figure(17)
plot(sig)
hold on
plot(outwf-sig,'r')
hold off
title(['Input and error signal of inverse filtered output with WF']);

% Listening
% sound([sig; sigout; sigoutnoise; siginv; outfreq; outwf],fs)






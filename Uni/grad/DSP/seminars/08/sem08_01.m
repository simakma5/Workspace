close all; clear; clc;
%% Input signal
load sink.mat;
clean = sum(sink(1:3,:));
clean = clean(:);  % force column vector

K = 0.3;
noise = loadbin("nc2.bin");
noise = K*noise(1:length(clean));

Ps = mean(clean.^2);
Pn = mean(noise.^2);
SNRin = 10*log10(Ps/Pn);
% % The same result can be obtained using a native Matlab function
% SNRin = snr(clean,noise);
% disp(SNRmatlab)

sig = clean + noise;

%% Output signal through Wiener filter
slen=length(sig);
wlen=512 ;
wstep=wlen/2 ;
wnum=(slen-wlen)/wstep+1 ;
w=hamming(wlen);

Ss = myPSD(clean,wlen,wstep);
Sn = myPSD(noise,wlen,wstep);
Hwf = Ss./(Ss+Sn);

% figure(99)
% subplot(311);
% plot(10*log10(Ss))
% axis tight
% title('Ss')
% 
% subplot(312)
% plot(10*log10(Sn))
% axis tight
% title('Sn')
% 
% subplot(313)
% plot(10*log10(Hwf))
% axis tight
% title('Hwf = Ss/(Ss+Sn)')

% The following is an adjusted OLA method for Wiener filtration
% Zero-initialization of output signal
out = zeros(slen,1);
out_dct = zeros(slen,1);

% Main cycle
for i=1:wnum

  % short-time frame selection
  ii=(i-1)*wstep+1;
  jj=(i-1)*wstep+wlen;

  frame=sig(ii:jj).*w;

  Sframe = fft(frame);
  Syframe = Hwf.*Sframe;
  yframe = ifft(Syframe);

  Sframe_dct = dctxc1(frame);
  Syframe_dct = Hwf.*Sframe;
  yframe_dct = ifft(Syframe);
  
  % Addition to output signal (implementation of OLA with general window)
  out(ii:jj)=out(ii:jj)+yframe;
  out_dct(ii:jj)=out_dct(ii:jj)+yframe_dct;

end

% Normalization of output signal
scale2 = sum(w(1:wstep:wlen));
out = out./scale2;
out_dct = out_dct./scale2;

% SNR of the output signal
Psw = mean(out.^2);
Pnw = mean((out-clean).^2);
SNRout = 10*log10(Psw/Pnw);

Psw_dct = mean(out_dct.^2);
Pnw_dct = mean((out_dct-clean).^2);
SNRout_dct = 10*log10(Psw_dct/Pnw_dct);

% Clean, input and output signals plotted
figure(1); clf;
subplot(411);
plot(clean)
title('The original signal without noise');

subplot(412);
plot(sig);
title(['Input signal, SNR: ' num2str(SNRin)]);

subplot(413);
plot(out) ;
title(['Output signal, SNR: ' num2str(SNRout)]);
% Note: the dampened edges are caused by the framing proccess
% Could be solved e.g., by padding the signal with zeros from both sides

subplot(414);
plot(out_dct) ;
title(['Output signal using DCT, SNR: ' num2str(SNRout_dct)]);

%% Noise dampening using amplitude spectrum difference
slen=length(sig);
wlen=512 ;
wstep=wlen/2 ;
wnum=(slen-wlen)/wstep+1 ;
w=hamming(wlen);

out = zeros(slen,1);
Nnoise = averageNoise(noise,wlen,wstep);

% Main cycle
for i=1:wnum

  % short-time frame selection
  ii=(i-1)*wstep+1;
  jj=(i-1)*wstep+wlen;

  sframe=sig(ii:jj).*w;
  Sframe = fft(sframe);
  Syframe = abs(Sframe)-Nnoise;

  for k=1:wlen
      % amplitude spectrum can't be negative -> fix if it occurs
      if Syframe(k) <= 0
          Syframe(k) = 0;
      % otherwise add phase
      else
          Syframe(k) = Syframe(k)*exp(1i*angle(Sframe(k)));
      end
  end
  yframe = ifft(Syframe);
  
  % Addition to output signal (implementation of OLA with general window)
  out(ii:jj)=out(ii:jj)+yframe;

end

% Normalization of output signal
scale2 = sum(w(1:wstep:wlen));
out = real(out./scale2);

% SNR of the output signal
Psw = mean(out.^2);
Pnw = mean((out-clean).^2);
SNRout = 10*log10(Psw/Pnw);

% Clean, input and output signals plotted
figure(2); clf;
subplot(311);
plot(clean)
title('The original signal without noise');

subplot(312);
plot(sig);
title(['Input signal, SNR: ' num2str(SNRin)]);

subplot(313);
plot(out) ;
title(['Output signal using amp. spec. diff., SNR: ' num2str(SNRout)]);
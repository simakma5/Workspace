%% AMPLITUDE SPECTRE SUBTRACTION PRINCIPLE TO FILTER OUT THE NOISE
% CV8 continuation here

% EXAMPLE OF OLA-GENERAL IMPLEMENTATION 
%  (without any modification)

%addpath H:\VYUKA\DSP\m
%addpath H:\VYUKA\DSP\signaly

% addpath /home/pollak/www/html/vyu/b2m31dsp/m
% addpath /home/pollak/www/html/vyu/b2m31dsp/signaly

% Loading of clean signal
%%%%%%%%%%%%%%%%%%%%%%%%%
%clean=loadbin('vm0.bin');  fs=16000 ;
%clean=clean(:);

load('sink.mat');

clean = sum(sink(1:3,:)); % clean signal suma prvnich 3 radky

sig = clean.'; %udela sloupcovy z cehokoliv celan = clean(:);

noise = loadbin('nc2.bin');
K = 0.3;
slength= length(sig);
noise = K*noise(1:slength,:);
sig = sig+noise;

Pcl = mean(clean.^2);
Pnc = mean(noise.^2);

SNR = 10*log10(Pcl/Pnc);

% Computation of short-time frame amount (50% overlap)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
slen=length(sig);
wlen=512 ;
wstep=wlen/2 ;
wnum=floor((slen-wlen)/wstep+1) ;

w=hamming(wlen);

% Vykresleni vstupniho signalu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf ;
subplot(211);
plot(sig);
title('Input signal');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zero-initialization of output signal
Noise_guess = zeros(wlen,wnum);

% NOISE CYCLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:wnum

  % short-time frame selection
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ii=(i-1)*wstep+1;
  jj=(i-1)*wstep+wlen;

  frame=noise(ii:jj).*w;
%we have to add the wiener filtration and spektralni vykonove hustoty;
%PSD pres pwelch or if pwelch if it does not sit; we need twosided spectrum
%from welch since from fft will be used twosided spectrum
% try to declare the welch
Sframe = fft(frame);
Sframe_abs = abs(Sframe);

  
  % No modification
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Addition to output signal (implementation of OLA with general window)
  Noise_guess(:,i)=Sframe_abs ; %poscitani framu

end

Noise_guess = mean(Noise_guess,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zero-initialization of output signal
out = zeros(slen,1);


% Main cycle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:wnum

  % short-time frame selection
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ii=(i-1)*wstep+1;
  jj=(i-1)*wstep+wlen;

  frame=sig(ii:jj).*w;
%we have to add the wiener filtration and spektralni vykonove hustoty;
%PSD pres pwelch or if pwelch if it does not sit; we need twosided spectrum
%from welch since from fft will be used twosided spectrum
% try to declare the welch
Sframe = fft(frame);
phase_frame = angle(Sframe);
Sframe_abs = abs(Sframe);
subtraction_fft = Sframe_abs-Noise_guess;
subtraction_fft = abs(subtraction_fft); %FW usmerneni
%subtraction_fft(subtraction_fft<0) = 0; %HW dvouceste usmereni
subtraction_fft = subtraction_fft.*exp(1j*phase_frame);




  
  % No modification
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  yframe = real(ifft(subtraction_fft));

  
  % Addition to output signal (implementation of OLA with general window)
  out(ii:jj)=out(ii:jj)+yframe ; %poscitani framu

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalization of output signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scale2 = sum ( w(1:wstep:wlen) ) ;
out = out ./ scale2 ;

subplot(212);
plot(out) ;
title('Output signal');

Pcl = mean(clean.^2);
out = out.';
yn = out(wlen:end-wlen)-clean(wlen:end-wlen);
Pyn = mean(yn.^2);

SNRy = 10*log10(Pcl/Pyn);

title(['Outout signal, SNR = ' num2str(SNRy,'%2.f') 'dB']);

SNRE = SNRy - SNR

%% Nonstationary signlas

% for nonstationary WF doesnt work but SPECTRAL SUBTRACTION WORKS FOR
% NONSTATIONARY SIGNAL (SINCE WE CAN NOT PREEMPTIVELY CALC WIENER)

% FOR NONSTATIONARY WE CAN CALCULATE ADAPTIVE WF, we CALCULATE Hwf for
% every separate frame

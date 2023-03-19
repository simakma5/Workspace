%DSP cheatsheet

close all; clc; clear all;

% signal
s = loadbin("SA001S01.CS0"); %sth random
fs = 16000;

%% Cepstrum (cv 6)
wlen = 512;
% wlen = floor(0.032 * fs); % 0.032 s
cp=12;
cepr = vrceps(s, 1, cp, wlen, wlen/2); %real cepstrum

%cepstral distance 
cepr1 = cepr;
cepr2 = cepr;
cdist = cde(cepr1, ceprf, cp); %%%%%%%%%%%%%%%%%%%%%%%%%nejde????
%(+compute mean )
meancd = mean(cdist);

%cepstr diff
ceprd = diff(cepr, 1, 1);
%distance
cd = sqrt(sum(ceprd.^2, 2)); %for changes in cepstrum !!


%% Filter (cv 6)
M=30;
Wc = 0.9; %LP cutoff frequency -- at 90% of freqs !!!!!!

b = fir1(M, Wc);
sf = filter(b, 1, s);

%% Spectrum (cv 6)
spectrogram(s, wlen, [],[], fs, 'yaxis')
pwelch(s, wlen, [],[], fs)

%% DCT (cv7)
dct1_a = dctxc1(s);
dct2_a = dctxc2(s);


%computed by dft
x_sym_1 = [x, x(x_len-1:-1:2)] %mirrored (excluding k = 0, N-1)
dct1_b = real(fft(x_sym_1));

x_sym_2 = [x, fliplr(x)] %mirrored (excluding k = 0, N-1)
dct2_b = real(fft(x_sym_2));

%% cepstrum from DFT and DCT (cv 7)
%rceps from DFT
dft1 = fft(s);
logdft1 = log(abs(dft1));
ceps3 = real(ifft(logdft1));
plot(ceps3, "DisplayName", "from DFT");

%rceps from DCT iDCT (DCT and iDCT same -- without constant)
dct1 = dctxc1(s);
logdct1 = log(dct1);
%ceps2 = real(idctxc1(logdct1(1:length(logdct1)/2))); %???? (ma to byt nejak z pulky vzorku, [ma byt pouzito DCT2 -- ok])
ceps2 = real(idctxc1(logdft1(1:length(logdft1)/2+1)));

%% SNR (cv 8)
u = 0.01 * randn(size(s))

P_s = mean(s.^2);
P_u = mean(u.^2);
snr = 10*log10(P_s/P_u)

%% his template (cv 8)

celan = s;
noise = u;
sig = s + u;

% Computation of short-time frame amount (50% overlap)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
slen=length(sig);
wlen=512 ;
wstep=wlen/2 ;
wnum=(slen-wlen)/wstep+1 ;

w=hamming(wlen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zero-initialization of output signal
out = zeros(slen,1);

% Main cycle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gs = [];
Gn  = [];

for i=1:wnum

  % short-time frame selection
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ii=(i-1)*wstep+1;
  jj=(i-1)*wstep+wlen;

  frame=sig(ii:jj).*w;
  
  % No modification 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --- DO STH WITH FRAME AND OUT
  yframe = frame ;
  

  Ss = abs(dctxc1(clean(ii:jj).*w)).^2 /wlen;
  Sn = abs(dctxc1(noise(ii:jj).*w)).^2 /wlen;

  %PWELCH of clean signal !!!!!!!!!!!!!!!
      Gs(:, i) = Ss;
      Gn(:, i) = Sn;
  
  %outframe = ifft(outframe);

  % Addition to output signal (implementation of OLA with general window)
%   out(ii:jj)=out(ii:jj)+yframe ;

end

%% WIENER (cv 8)
%neorig with sqrt() -- lower 
H_wf = sqrt( Ss ./ (Ss + Sn) ); %FOR NOISE FILTERING (after cycle)

%% SPECTRAL DIFF (cv 8)
% in that cycle !!!

  Sx = fft(frame);
  phases = angle(Sx);

   Y = abs(Sx) - sqrt(wlen * Sn); 

   %one-way rectifying !! (stronger, but more distortion)
%    Y(Y<0) = 0;
   %two-way rectifying !! (weaker but, less distortion)
   Y(Y<0) = -Y(Y<0);

   y_this_segment = real(ifft(Y .* exp(1i * phases)));

   

   %% only deconv (cv 10)
   Y_fr_inv_filtred = (1 ./ H) .* Y_fr;

    %% WIENER WITH DECONVOLUTION (cv 10)
    %fisrt
    S_uu = S_uu  +  abs(fft(frame)).^2 /wlen   / wnum; %G_noise

    %then
   Y_fr = fft(frame);
   S_yy = abs(fft(frame)).^2 /wlen;

   H_wf = (1 ./ H) .* ( abs(H).^2 ./ (abs(H).^2 + (S_uu ./ S_ss) ) ); %%%% WIENER !!
   Y_fr_inv_filtred = H_wf .* Y_fr; 



    %% % Analýza vlastních komponent signálu a KLT (cv 11)

    % Odhadněte kovarianční matici z 20 realizací daného signálu (viz níže, pro jednotlivé realizace uvažujte segmentaci oknem délky 200 vzorků s krokem 50 vzorků)
wstep = 50;
wlen = 200;

% generate our signals as segments of sig
SS = [];
for i = 1:20
    ii = (i-1)*wstep+1;
    jj = (i-1)*wstep+1 + wlen;

    SS = [SS; sig(ii:jj)];

end

Rss = SS' * SS; %pres realizace ???? -- 'mean' pres ty k. cleny nasich segmentu --> 200x200


% eigenvects and eigenvals
[V, D  ] = eig(Rss);
%greatest first (flip only)
d = diag(D); 
d = flipud(d);
V = fliplr(V);

d = diag(D);

% KLT
x = SS(1,:)'; %vec

Wklt = V'; %RADKY TRANSFORMACNI MATICE JSOU VL VEKTORY !!

%transformovani
Xklt = Wklt*x;


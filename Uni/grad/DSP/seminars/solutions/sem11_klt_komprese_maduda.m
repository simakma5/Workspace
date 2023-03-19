clc; 
% clear all; close all;

% Setup of directory with signals for different OS
if exist('/home/pollak/vyuka/DSP/m'),  % it is the Linux OS
    addpath /home/pollak/vyuka/DSP/m
end
if exist('/home/pollak/vyuka/DSP/signaly'),  % it is the Linux OS
    addpath /home/pollak/vyuka/DSP/signaly
end

if exist('H:\VYUKA\DSP\m'),
    addpath H:\VYUKA\DSP\m
end
if exist('H:\VYUKA\DSP\signaly'),
    addpath H:\VYUKA\DSP\signaly
end

if exist('C:\usr\pollak\vyuka\DSP\m')
    addpath C:\usr\pollak\vyuka\DSP\m
end
if exist('C:\usr\pollak\vyuka\DSP\signaly')
    addpath C:\usr\pollak\vyuka\DSP\signaly
end



%% s1 + noise
fs=8000;
f1=821; w1=821*2*pi;
slen=2000;
t=linspace(0,((slen-1)/fs),slen);
K=0.01;
s1=sin(w1*t)+K*randn(1,slen);
s1=s1-mean(s1);


%% s2 + noise

f21=532; w21=f21*2*pi; 
f22=640; w22=f22*2*pi;
s2=sin(w21*t)+sin(w22*t)+K*randn(1,slen);
s2=s2-mean(s2); %centering


%% r1
r=loadbin('vm0.bin');
slen=length(r);
fsr1=16000;
r1=r+K*randn(slen,1);
r1=r1-mean(r1); %centering

%% Selection of signal

% s=s1;    %%% one sinusoid with white noise
% s=s2;    %%% two sinusoids with white noise
s=r1;    %%% speech sginal with white noise


%%segmentace

slen=length(s);
wlen=200;
wstep=50;
wnum=(slen-wlen)/wstep;
RR=zeros(wnum, wlen); i=1;


for start=0:wstep:wnum*(wstep+1)
      RR(i,:)=s(start+1:start+wlen);
      i=i+1;
end;

%transponované datové matice
RRt=RR.';

%kovariační matice
COV1=RRt*RR;

%vlastní čísla a vektory
[eVect_s, eVal_s]=eig(COV1); 
eVect_s=flip(eVect_s,2);
eVal_s=flip(eVal_s,2);
eigVal_s=reshape(eVal_s,[],1);
eigVal_s=eigVal_s(eigVal_s~=0);




%% Porovnání kompresních vlastností
%KLT

Wklt = eVect_s';
Wkltinv = Wklt';

sig=RR(10,:);
Klt=Wklt*sig';

%DFT
Dft=abs(fft(sig)); %twosided
Dft_complex=fft(sig);
Phase=angle(Dft_complex);

% Creation of onesided DFT spectrum
Dft_single=Dft(1:length(Dft)/2);
Dft_single(2:end)=2*Dft_single(2:end); %singlesided

%DCT
Dct=dct(sig);

%Kumulativní součty koeficientů 
csum_dft=cumsum(Dft_single); dft95=0.95*max(csum_dft);
csum_klt=cumsum(abs(Klt)); klt95=0.95*max(csum_klt);
csum_dct=cumsum(abs(Dct)); dct95=0.95*max(csum_dct);
%Počet komponent pro rekonstrukci signálu s 95% přesností
component_nr_klt=find(csum_klt>klt95,1);
component_nr_dft=find(csum_dft>dft95,1);
component_nr_dct=find(csum_dct>dct95,1);


figure(1); subplot 321; stem(Klt); hold on; xline(component_nr_klt,'r','Linewidth',1.0); title('KLT obraz'); hold off
subplot 323; stem(Dft_single); hold on; xline(component_nr_dft,'r','Linewidth',1.0); title('DFT spektrum'); hold off
subplot 325; stem(Dct); hold on; xline(component_nr_dct,'r','Linewidth',1.0); title('DCT spektrum'); hold off

subplot 322; stem(csum_klt/max(csum_klt)); hold on; yline(0.95); xline(component_nr_klt,'r','Linewidth',1.0); title('KLT obraz'); title('KLT kumulativní součet'); hold off
subplot 324; stem(csum_dft/max(csum_dft)); hold on; yline(0.95); xline(component_nr_dft,'r','Linewidth',1.0); title('DFT spektrum'); title('DFT kumulativní součet'); hold off
subplot 326; stem(csum_dct/max(csum_dct)); hold on; yline(0.95); xline(component_nr_dct,'r','Linewidth',1.0); title('DCT spektrum'); title('DCT kumulativní součet'); hold off

%% zpětná rekonstrukce
iDct=zeros(size(Dct)); iDct(1:component_nr_dct)=Dct(1:component_nr_dct);
iDft=zeros(size(Dft)); iDft(1:component_nr_dft)=Dft(1:component_nr_dft); iDft(end-component_nr_dft:end)=Dft(end-component_nr_dft:end);
iKlt=zeros(size(Klt)); iKlt(1:component_nr_klt)=Klt(1:component_nr_klt);

backsig_dft=real(ifft(iDft.*exp(j*Phase)));
backsig_dct=idct(iDct);
backsig_klt=Wkltinv*iKlt;

figure(2);
subplot 411; plot(sig); title('Původní signál');
subplot 412; plot(backsig_klt); title(['KLT rekonstrukce; počet komponent: ' num2str(component_nr_klt)]);
subplot 413; plot(backsig_dft); title(['DFT rekonstrukce; počet komponent: ' num2str(2*component_nr_dft)]);
subplot 414; plot(backsig_dct); title(['DCT rekonstrukce; počet komponent: ' num2str(component_nr_dct)]);
component_nr_dct
component_nr_dft=2*component_nr_dft
component_nr_klt




    

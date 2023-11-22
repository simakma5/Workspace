clear variables
close all

Tp=25e-6;   % pulse length

% load('DataCv7b.mat')
load('DataCv7a.mat')

%Processing
  % pulse compression
sSumComp=filter(flipud(conj(s0.*hamming(length(s0))))',1,sSum,[],1);
  % MTI processing
sSumCompMTI=filter([1 -2 1]/4,1,sSumComp,[],2);
  % MTD processing
NMTD=32;

% Two types of MTD filters - uncomment to compare
% example of non-windowed MTD filter (rect)
% hMTD=repmat(ones(NMTD,1)',NMTD,1).*exp(1i*2*pi*(0:NMTD-1)'*(0:NMTD-1)/((NMTD)));
% example of Hamming windowed MTD filter
hMTD=repmat(hamming(NMTD)',NMTD,1).*exp(1i*2*pi*(0:NMTD-1)'*(0:NMTD-1)/((NMTD)));

for ii=1:NMTD
  sSumCompMTD{ii}=filter(hMTD(ii,:),1,sSumComp,[],2); %#ok<SAGROW>
end


% Display ampl input/kompress/MTI filterred
figure('name','Signal2D')
subplot(3,1,1)
imagesc(20*log10(abs(sSum)),[0,100])
colormap('jet')
hb=colorbar();
ylabel(hb,'A [dBLSB]')
grid on
title ('Input signal')
xlabel('Azim quanta')
ylabel('Range quanta')
subplot(3,1,2)
imagesc(20*log10(abs(sSumComp)),[0,100])
colormap('jet')
hb=colorbar();
ylabel(hb,'A [dBLSB]')
grid on
title ('Compressed signal')
xlabel('Azim quanta')
ylabel('Range quanta')
subplot(3,1,3)
% imagesc(20*log10(abs(filter([1 -2 1],1,filter(flipud(conj(s0.*hamming(length(s0)))),1,sSum,[],1),[],2))))
imagesc(20*log10(abs(sSumCompMTI)),[0,100])
colormap('jet')
hb=colorbar();
ylabel(hb,'A [dBLSB]')
grid on
title ('MTI filtered signal')
xlabel('Azim quanta')
ylabel('Range quanta')
%%
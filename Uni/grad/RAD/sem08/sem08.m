close all; clear; clc; addpath(genpath(fullfile([pwd '\Uni' '\grad' '\RAD'])))

%% Initialization
load DataCv7a.mat
% load DataCv7b.mat
Tp = 25e-6;

%% Processing (initial part is the same processing as in sem07)
% Pulse compression (hamming windowing + matched filtering)
sSumComp = filter(flipud(conj(s0.*hamming(length(s0))))', 1, sSum, [], 1);

% MTI processing (double suppresion / three pulse canceller)
sSumCompMTI = filter([1 -2 1]/4, 1, sSumComp, [], 2);

% MTD processing (two types -- uncomment to compare)
NMTD = 128;    % number of pulses processed
% non-windowed (rectangle) vs. Hamming-windowed MTD filter
% hMTD = repmat(ones(NMTD,1)', NMTD, 1).*exp(1i*2*pi*(0:NMTD-1)'*(0:NMTD-1)/NMTD);
hMTD = repmat(hamming(NMTD)', NMTD, 1).*exp(1i*2*pi*(0:NMTD-1)'*(0:NMTD-1)/NMTD);
sSumCompMTD = cell(1, NMTD);
for pulse = 1:NMTD
    sSumCompMTD{pulse} = filter(hMTD(pulse,:), 1, sSumComp, [], 2);
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


figure(100)
mesh(20*log10(abs(sSumCompMTI(:,3:end))))
xlabel('pulse number(azim. quanta)')
ylabel('range quanta [-]')
title('MTI - range vs. pulse index')
%%
sSumCompMTIacc=sum(sSumCompMTI(:,3:end),2);
figure(20)
plot(20*log10(abs(sSumCompMTIacc)))
xlabel('range quanta [-]')
title('MTI coherent integrated') 
sSumCompMTIaccncoh=sum(abs(sSumCompMTI(:,4:end)),2);
figure(21)
plot(20*log10(sSumCompMTIaccncoh))
xlabel('range quanta [-]')
title('MTI non coherent integrated')
%CFAR
CFARinterval=10;
guardinterval=1;
detection_vector=zeros(1,length(sSumCompMTIaccncoh)-CFARinterval-1);
for ii=CFARinterval+1:length(sSumCompMTIaccncoh)-CFARinterval-1
    background_estimate_1=sum(sSumCompMTIaccncoh(ii-CFARinterval:ii-1-guardinterval));
    background_estimate_2=sum(sSumCompMTIaccncoh(ii+guardinterval+1:ii+CFARinterval));
    background_estimate = max([background_estimate_1,background_estimate_1]);
    if(sSumCompMTIaccncoh(ii) > background_estimate)
        disp(['Target detected at ', num2str(ii), 'sample'])
        detection_vector(ii)=1;
    end
end
figure(22)
plot(detection_vector)
xlabel('range quanta [-]')
title('MTI non coherent integrated - target threshold')

target_index=[];
detect_target=0;
accum_indexes=0;
count_detections=0;
for kk=1:length(detection_vector)
    if(detection_vector(kk)==1)
        detect_target=1;
        accum_indexes=accum_indexes+kk;
        count_detections=count_detections+1;
    end
    if(detection_vector(kk)==0)&&(detect_target==1)
        target_index=[target_index,accum_indexes/count_detections];
        detect_target=0;
        accum_indexes=0;
        count_detections=0;
    end
end
target_distances=target_index*3e8/fs/2+Rmin-Tp*3e8/2;     
disp(['estimated target distances: ',num2str(target_distances), ' m'])        
        
%%
sSumCompMTDacccoh = zeros(length(sSum),length(sSumCompMTD));
for ii = 1:NMTD
    sSumCompMTDacccoh(:,ii) = sum(sSumCompMTD{1,ii},2);
end
figure(30)
mesh(20*log10(abs(sSumCompMTDacccoh)))
xlabel('Doppler bank frequency index')
ylabel('range quanta [-]')
title('MTD - coherent integration')

sSumCompMTDaccncoh = zeros(length(sSum),length(sSumCompMTD));
for ii = 1:NMTD
    sSumCompMTDaccncoh(:,ii) = sum(abs(sSumCompMTD{1,ii}),2);
end
figure(31)
mesh(20*log10((sSumCompMTDaccncoh)))
xlabel('Doppler bank frequency index')
ylabel('range quanta [-]')
title('MTD - noncoherent integration')

%%
for ii=1:NMTD
  sSumCompMTIMTD{ii}=filter(hMTD(ii,:),1,sSumCompMTI,[],2); %#ok<SAGROW>
end
sSumCompMTIMTDaccncoh = zeros(length(sSum),length(sSumCompMTD));
for ii = 1:NMTD
    sSumCompMTIMTDaccncoh(:,ii) = sum(abs(sSumCompMTIMTD{1,ii}),2);
end
figure(40)
xscale = [(0:NMTD-1)/NMTD*PRF];
% xscale = xscale*3e8/fc/2*3.6;
yscale = Rmin-Tp*3e8/2+[(0:length(sSumCompMTIMTDaccncoh)-1)*3e8/fs/2];
% mesh(xscale,yscale,20*log10((sSumCompMTIMTDaccncoh)))
surf(xscale,yscale,20*log10((sSumCompMTIMTDaccncoh)))

xlabel('Doppler frequency [Hz]')
% xlabel('radial speed [kmph]')
% xlabel([num2str((0:NMTD-1)/NMTD*PRF)])
ylabel('Slant range [m]')
title('MTI+MTD - noncoherent integration')

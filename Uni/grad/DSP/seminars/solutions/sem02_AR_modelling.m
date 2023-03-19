% LPC - AR modelling
%%%%%%%%%%%%%%%%%%%%%%%%

% 1) Analyzed signal

frame = loadbin ( 'vm0-512.bin' ) ;

% Computation of AR model parameters

% The choice of the order
% p=6;
% p=10;
p=16;
[a, Ep]=lpc(frame,p);


% 1) Filtering by analyzing filter - obtaining the error signal

figure(1) ;

subplot(211) ;
plot(frame) ;
title('Original signal')

en=filter(a,1,frame) ;

subplot(212)
plot(en)
title('Error signal') ;


% 2) Modelling with the excitation by error signal

figure(2) ;

subplot(211) ;
plot(en) ;
title('Excitation of AR model - error signal')

y1=filter(1,a,en) ;

subplot(212)
plot(y1,'r')
hold on
plot(frame,'--')
hold off
title('Modelled and original signal') ;

% 3) Modelling with the excitation by artificailly generated pulses
%    BAD POWER ADJUSTMENT

figure(3) ;

To=120 ;
un=[1 zeros(1,To-1) 1 zeros(1,To-1) 1 zeros(1,To-1) 1 zeros(1,To-1) 1 zeros(1,To-1)] ;


subplot(211) ;
plot(un) ;
title('Excitation of AR model - artificial pulses - BAD power')

y2=filter(1,a,un) ;

subplot(212)
plot(y2,'r')
hold on
plot(frame,'--')
hold off
title('Modelled and original signal') ;


% 4) Modelling with the excitation by artificailly generated pulses
%    CORRECT SETTING OF THE POWER

figure(4) ;

Pun = mean(un.^2) ;
un2=un./sqrt(Pun) ;


subplot(211) ;
plot(un2) ;
title('Excitation of AR model - artificial pulses - CORRECT power')

y3=filter(sqrt(Ep),a,un2) ;

subplot(212)
plot(y3,'r')
hold on
plot(frame,'--')
hold off
title('Modelled and original signal') ;




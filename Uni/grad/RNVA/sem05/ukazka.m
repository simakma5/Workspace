% s = [ones(1,5) zeros(1,15)]; %obdelnik
sp = [1 -1 -1 1 -1 1 1]; %PN posloupnost m-sequence m=3 [3 2 0]
% sp = [1 -1 -1 1 1 1 -1]; %PN posloupnost m-sequence m=3 [3 1 0]
% sp = [ 1 1 1 1 1 -1 1 1 -1 -1 1 1 1 -1 -1 -1 -1 1 1 -1 1 -1 1 -1 -1 1 -1 -1 -1 1 -1];
% sp = [-1 -1 1 -1 -1 -1 -1 -1 -1 1 1 -1 -1 1 1 1 -1 1 1 1 1 -1 -1 1 1 -1 -1 1 1 1 1];

sampling_rate=10;

Barker_codes % nacteni Barkerovych kodu (9 variant) ze souboru Barker_codes.m

%% Selected Barker code wavwform and crosscorrelation function
s=[0, Barker11, 0];
N = length(s);
R = xcorr(s,s);
figure(11)
plot(-N+1:N-1,R)

% sig=[];
% for i=1:length(s)
%     sig=[sig, s(i)*ones(1,sampling_rate)];
% end
sig=repmat(s,sampling_rate,1);
sig=sig(:)';
figure(1)
plot(sig)

%% Overview of all Barker codes
figure(12)
subplot(3,3,1)
s=[0, Barker2_1, 0];
N = length(s);
R = xcorr(s,s);
plot(-N+1:N-1,R)
title('Barker 2 var.1')
figure(2)
subplot(3,3,1)
sig=repmat(s,sampling_rate,1);
sig=sig(:)';
plot(sig)
title('Barker 2 var.1')


figure(12)
subplot(3,3,2)
s=[0, Barker2_2, 0];
N = length(s);
R = xcorr(s,s);
plot(-N+1:N-1,R)
title('Barker 2 var.2')
figure(2)
subplot(3,3,2)
sig=repmat(s,sampling_rate,1);
sig=sig(:)';
plot(sig)
title('Barker 2 var.2')


figure(12)
subplot(3,3,3)
s=[0, Barker3, 0];
N = length(s);
R = xcorr(s,s);
plot(-N+1:N-1,R)
title('Barker 3')
figure(2)
subplot(3,3,3)
sig=repmat(s,sampling_rate,1);
sig=sig(:)';
plot(sig)
title('Barker 3')

figure(12)
subplot(3,3,4)
s=[0, Barker4_1, 0];
N = length(s);
R = xcorr(s,s);
plot(-N+1:N-1,R)
title('Barker 4 var.1')
figure(2)
subplot(3,3,4)
sig=repmat(s,sampling_rate,1);
sig=sig(:)';
plot(sig)
title('Barker 4 var.1')

figure(12)
subplot(3,3,5)
s=[0, Barker4_2, 0];
N = length(s);
R = xcorr(s,s);
plot(-N+1:N-1,R)
title('Barker 4 var.2')
figure(2)
subplot(3,3,5)
sig=repmat(s,sampling_rate,1);
sig=sig(:)';
plot(sig)
title('Barker 4 var.2')

figure(12)
subplot(3,3,6)
s=[0, Barker5, 0];
N = length(s);
R = xcorr(s,s);
plot(-N+1:N-1,R)
title('Barker 5')
figure(2)
subplot(3,3,6)
sig=repmat(s,sampling_rate,1);
sig=sig(:)';
plot(sig)
title('Barker 5')

figure(12)
subplot(3,3,7)
s=[0, Barker7, 0];
N = length(s);
R = xcorr(s,s);
plot(-N+1:N-1,R)
title('Barker 7')
figure(2)
subplot(3,3,7)
sig=repmat(s,sampling_rate,1);
sig=sig(:)';
plot(sig)
title('Barker 7')

figure(12)
subplot(3,3,8)
s=[0, Barker11, 0];
N = length(s);
R = xcorr(s,s);
plot(-N+1:N-1,R)
title('Barker 11')
figure(2)
subplot(3,3,8)
sig=repmat(s,sampling_rate,1);
sig=sig(:)';
plot(sig)
title('Barker 11')

figure(12)
subplot(3,3,9)
s=[0, Barker13, 0];
N = length(s);
R = xcorr(s,s);
plot(-N+1:N-1,R)
title('Barker 13')
figure(2)
subplot(3,3,9)
sig=repmat(s,sampling_rate,1);
sig=sig(:)';
plot(sig)
title('Barker 13')

%% One period of PRN code
s = sp;
s = [sp, sp, sp, sp,sp];
Np = length(s);
R2 = xcorr(s,s);
figure(3)
plot(-Np+1:Np-1,R2)
title('xcorr PRN')
%%
Rp = korelace(s,s);
figure(4)
plot(0:length(s)-1,Rp)
title('cyklicka korelace PRN') 

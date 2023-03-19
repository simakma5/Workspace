clear; clc; close all;

%%
% A = [1 0.5; 0.5 1];
% [V, D] = eig(A);
% disp(V);disp(diag(D)');

%%
f = 821;
f1 = 532;
f2 = 640;
fs = 8e3;
slen = 2e3;
scale = 1e-3;
t = (0:slen-1)/fs;

s1 = sin(2*pi*f*t) + scale*randn(1,slen);
s2 = sin(2*pi*f1*t) + sin(2*pi*f2*t) + scale*randn(1,slen);
r1 = loadbin("vm0.bin")';
r1 = r1 + scale*randn(size(r1));

% sig = s1; fs = 8e3;
sig = s2; fs = 8e3;
% sig = r1; fs = 16e3;

wlen = 200;
wstep = 50;
SS = [];
for i = 1:20
    ii = (i-1)*wstep+1;
    jj = (i-1)*wstep+wlen;

    SS = [SS; sig(ii:jj)];  % strednik za SS udava spravne formatovani
end

Rss = SS'*SS;               % Rss = SS^{\mathrm T} \cdot SS

%%
[V,D] = eig(Rss);
d = diag(D);                % vytahnu jenom eigenvalues
d = flipud(d);              % dam na zacatek nejvyznamnejsi hodnoty
V = fliplr(V);              % dam na zacatek nejvyznamnejsi hodnoty

figure(1)
subplot(211);
plot(SS(1,:))

subplot(212);
stem(d)

figure(2)
for k = 1:4
    subplot(['41' num2str(k)]);
    plot(V(:,k))
end

%%
x = SS(1,:)';
% x = r1(1:wlen)';                % wrong transformation basis demonstration

Wklt = V';
Xklt = Wklt*x;
Xdft = fft(x);
Xdct = dctxc1(x);

figure(3)
subplot(311);
plot(abs(Xklt),'*-')             % abs to get real values

subplot(312);
plot(abs(Xdft),'*-')             % already real - regular abs

subplot(313);
plot(abs(Xdct),'*-')             % already real - regular abs

%%
Xdft = Xdft(1:end/2+1);
Xdft(2:end-1) = 2*Xdft(2:end-1);

% vypocet kumulativnich sum
csXklt = cumsum(abs(Xklt))/sum(abs(Xklt));
csXdft = cumsum(abs(Xdft))/sum(abs(Xdft));
csXdct = cumsum(abs(Xdct))/sum(abs(Xdct));

figure(4)
subplot(321);
plot(abs(Xklt),'.-')

subplot(323);
plot(abs(Xdft),'.-')

subplot(325);
plot(abs(Xdct),'.-')

subplot(322);
plot(abs(csXklt), '.-')

subplot(324);
plot(abs(csXdft), '.-')

subplot(326);
plot(abs(csXdct), '.-')
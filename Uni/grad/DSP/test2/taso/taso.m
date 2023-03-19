%% DSP test 2
%% 1
sig = loadbin('frame-024.bin');
sig = sig.*hamming(length(sig));
Csig = rceps(sig);
plot(Csig(1:20))

%% 2
load sigs_2chan_01.mat;
samp = 16*8;
w = hamming(samp);

msc_val = mscohere(sig1, sig2, w, samp/2, samp, 16e3);
msc_mean = mean(msc_val)

%% 3
sig1 = loadbin('SA013S01.CSX');
sig2 = loadbin('SA013S01.CS0');

p = 16;
cp = 20;

ceps1 = vaceps(sig1, 1, p, cp, 512);
ceps2 = vaceps(sig2, 1, p, cp, 512);

d = cde(ceps1(:, 2:end), ceps2(:, 2:end), cp);
dist = mean(d)

%% 4
sig1 = loadbin('frame-001.bin');
sig2 = loadbin('frame-009.bin');

ceps1 = vrceps(sig1, 1, 12, 512);
ceps2 = vrceps(sig2, 1, 12, 512);

cp = 12;
d = cde(ceps1, ceps2, cp);
dist = mean(d)

%% 5
load sig_xy_04.mat
wlen = 512;
w = hamming(wlen);
Nfft = wlen;
[S , ff] = cpsd(x, y, w, wlen/2, fs, 'onesided');
plot(ff, unwrap(angle(S)))

%% 6
snoi = loadbin('SX010S01.CS0');
scl = loadbin('SA010S01.CS0');

SNR_db = snr(scl, snoi-scl)

%% 7
sig = loadbin('frame-011.bin');
Sdct = dct(sig);
pwrdct = sum(Sdct(1:55).^2);
pwrsig = sum(sig.^2);
pct = pwrdct / pwrsig * 100

%% 8
sig = loadbin('frame-003.bin');
sig = sig.*hamming(length(sig));
Sdct = dctxc2(sig);
Sdct(1:8)

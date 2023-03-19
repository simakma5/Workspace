fs = 16e3;

sig1 = loadbin('SA001S01.CS0');
slen = length(sig1);
sig2 = sig1(20+1:20+slen);

wlen = 1024;
wstep = wlen/4;
noverlap = wlen - wstep;
wnum = floor((slen - wlen)/wstep) + 1;

wlen2 = 256;
noverlap2 = 3*wlen2/4;
wstep2 = wlen2 - noverlap2;
wnum2 = floor((wlen - wlen2)/wstep2) + 1;
nfft = max([wlen2 256]);

MSCmtx = zeros(nfft/2+1,wnum);
% tic
for k = 1:wnum
    ii = (k-1)*wstep + 1;
    jj = (k-1)*wstep + wlen;
    frame1 = sig1(ii:jj);
    frame2 = sig2(ii:jj);

    [msc,ff] = mscohere(fram1,frame2,wlen2,noverlap2,nfft);
    MSCmtx(:,k) = msc;
end
% toc


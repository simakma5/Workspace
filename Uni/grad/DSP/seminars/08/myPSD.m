function S = myPSD(sig,wlen,wstep)

sig = sig(:);

slen=length(sig);
wnum=floor((slen-wlen)/wstep+1);

w=hamming(wlen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zero-initialization of output signal
X = zeros(wlen,wnum);

% Main cycle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:wnum

  % short-time frame selection
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ii=(i-1)*wstep+1;
  jj=(i-1)*wstep+wlen;

  frame=sig(ii:jj).*w;

  X(:,i) = abs(fft(frame)).^2/wlen;  % zapis do i-teho sloupce

end

S = mean(X,2);
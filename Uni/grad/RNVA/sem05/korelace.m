function Rp = korelace(s1,s2)
N = length(s1);

% varianta 1
% Rp = zeros(1,N);
% for t = 1:N
%     for k = 1:N
%         Rp(t) = Rp(t) + s1(mod(k+t-2,N)+1)*s2(k)/N;
%     end
% end

% varianta 2
% for t = 1:N
%     Rp(t) = sum(circshift(s1,[0 -t+1]).*s2)/N;
% end

% varianta 3
% reseni korelace jako konvoluce f(x) s conj(f(-x))
% Rp = ifft(fft(s1).*conj(fft(s2)))/N;

% varianta 4
Rp2 = xcorr([s1 s1], s2)/N; Rp = Rp2(2*N:3*N-1);
%


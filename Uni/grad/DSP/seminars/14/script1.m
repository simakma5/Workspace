sig = loadbin("nc1.bin");

x = sig + 0.01*randn(size(sig));
% z = qmfb();

N = 3;
stop = 10000/2^N;
y = zeros(length(x));
for i = 1:length(x)
    the = mean(abs(x(1:stop,i)));
    y(:,i) = x(:,i).*(abs(x(:,1)) > thr(i));
end
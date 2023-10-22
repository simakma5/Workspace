close all; clear; clc

%% Mandelbrot set
N = 500;
Niter = 500;
x = repmat(linspace(-2, 1, N), N, 1);
y = repmat(linspace(-2, 1.5, N)', 1, N);
z0 = x + 1i*y;
mandelbrot = ones(N);
z = z0;
for n = 1:Niter
    z = z.*z + z0;
    for i = 1:N
        for j = 1:N
            if abs(z(i,j)) < 2
                mandelbrot(i,j) = mandelbrot(i,j) + 1;
            end
        end
    end
end
mandelbrot = log(mandelbrot);

% Visualization
figure("Name", "Mandelbrot set visualization")
imagesc(mandelbrot)
xlabel("x domain")
ylabel("y domain")
title(['Visualization of the Mandelbrot set for N = ' num2str(N) ', Niter = ' num2str(Niter)])
colormap(flipud(hot))
close all; clear; clc

%%
% Calculation of golden ratio
N = 16;
nRange = 1:N-1;
goldenRatio = zeros(1,N-1);
for n = nRange
    goldenRatio(n) = myFibonacci(n+1)/myFibonacci(n);
end

% Interpolation using spline, pchip or makima (can be interchanged)
fineRange = 1:1e-2:N-1;
interpolation = pchip(nRange, goldenRatio, fineRange);

% Plot results
figure("Name", "Golden ratio")
plot(nRange, goldenRatio, 'o', 'Color', 'k')
hold on
plot(fineRange, interpolation, '-', 'Color', 'k')
yline((1+sqrt(5))/2, '--')
axis tight
grid on
grid minor
xlabel("N")
ylabel({'$\varphi$'},'Interpreter','latex')
legend("Golden ratio", "Spline interpolation")
title("Golden ratio as the ratio of two consecutive Fibonacci numbers")
hold off

% Experimental Data Analysis: Lecture 6
close all, clear all, clc

%% MATLAB EXAMPLE 1 %%
% testing the alpha error

n = 100; % number of iterations
alpha = 0.05; % set alpha level

X = rand(100,n); % generate random numbers

for i = 1:n,
    % independent t-test of randomly split samples into two groups
    [h, p(i)] = ttest2(X(1:50,i),X(51:100,i)); 
    
end

p(p < alpha) % the fake results
100*length(p(p < alpha))/n % percentage number of fake results
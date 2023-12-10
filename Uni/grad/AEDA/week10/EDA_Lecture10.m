% Experimental Data Analysis: Lecture 10
close all, clear all, clc

%% MATLAB EXAMPLE 1 %%
%  Area under curve

ES = 1; % provide effect size for 2 random distributions

% generate 2 random distributions
x1 = randn(1,1000);
x2 = randn(1,1000)+ES;
x = [x1 x2];

% generate group lables
y(1:1000)=0;
y(1001:2000)=1;
y = logical(y)';

% perform a logistic regression
modelparams = glmfit(x,y,'binomial','link','logit');

% compute the class assignment (probabilities) of the trained model for each data point
modelfitLR = glmval(modelparams,x,'logit');

% compute Receiver Operating Characteristic (ROC) curve and Area Under
% Curve (AUC)
[X,Y,T,AUC] = perfcurve(y,modelfitLR,'true');

% plot ROC curve
plot(X,Y, 'b','LineWidth',2)
set(gca, 'Box','Off')
ylabel('SEN'); 
xlabel('1-SPC');


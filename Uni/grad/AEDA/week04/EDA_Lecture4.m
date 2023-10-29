% Experimental Data Analysis: Lecture 4
close all, clear all, clc

%% MATLAB EXAMPLE 1 %%
% estimate correlations between maximum phonation time until voice break 
% and motor clinical score of speakers with Huntington's disease using both
% Pearson and Spearman correlation

% load data for  maximum phonation time until voice break
mptvb=[2.7650 3.6125 1.8325 15.2625 4.6025 4.5750 3.8575 7.2125 4.5350 1.1900 ...
      2.2075 0.3650 0.7875 1.1075 2.5000 10.4550 7.4325 1.1650 9.7525 1.4925 ...
      7.3825 2.0875 4.8775 6.3200 22.6225 11.2425 9.0325 1.3550 4.9825 ...
      3.2775 7.9975 13.5275];
% load data for motor clinical score
uhdrs = [21 11 12 7 5 11 11 12 11 17 6 26 18 16 10 3 3 11 7 7 5 8 2 5 8 0 2 7 ...
         10 4 11 3];
% make a figure
figure(),
plot(mptvb, uhdrs,'x','MarkerSize',10,'LineWidth',2)
ylabel ('UHDRS','FontSize',12, 'FontName','Georgia')
xlabel ('mptvb (-)','FontSize',12, 'FontName','Georgia')
set(gca, 'Box', 'off',  'FontName','Georgia')

% Pearson parametric correlation
[r,p] = corr(mptvb(:), uhdrs(:), 'type', 'Pearson') 

% Spearman non-parametric correlation
[r,p] = corr(mptvb(:), uhdrs(:), 'type', 'Spearman') % Spearman correlation

%% MATLAB EXAMPLE 2 %%
% test the normality of data based upon 4 measures extracted from
% sustained phonation of 34 healthy speakers

load data_phonation; % load data for 
MPT = phonation(:,1); % Maximum phonation time
F0SD = phonation(:,2); % Pitch variations
Shimmer = phonation(:,3); % Amplitude perturbation
DUV = phonation(:,4); % Degree of unvoiced segments

% Kolmogorov-Smirnov test
variable = F0SD;

% KS test including z-score normalization
[h,p]=kstest((variable-mean(variable))/std(variable)) 

% KS test through creation of cumulative distribution function
normdata=normcdf(variable,mean(variable),std(variable)); 
[h,p,ksstat] = kstest(variable,[variable,normdata],0.05)

% plotting the data
figure(), subplot(121)
[N,X] = hist(variable, 20); 
bar(X, N, 'b');

% plotting the CDFs
subplot(122)
[f,x_values] = ecdf(variable);
EmpCDF = plot(x_values,f);
hold on;
StandCDF = plot(sort(variable),sort(normdata),'r');
set(EmpCDF,'LineWidth',2);
set(StandCDF,'LineWidth',2);
legend([EmpCDF  StandCDF],'Empirical CDF','Standard Normal CDF','Location','SE');

%% MATLAB EXAMPLE 3 %%
% Tests of normality available in Matlab

normdata = randn(1,100); % generate random data for normal distribution
abnormdata = (randn(1,100)).^2 + 1;  % generate random data for normal distribution

% Chi-Square test
[h,p] = chi2gof(normdata,'cdf',@normcdf)
[h,p] = chi2gof(abnormdata,'cdf',@normcdf)

% Liliefors test
[h,p] = lillietest(normdata)
[h,p] = lillietest(abnormdata)

% Anderson-Darling test
[h,p] = adtest(normdata)
[h,p] = adtest(abnormdata)


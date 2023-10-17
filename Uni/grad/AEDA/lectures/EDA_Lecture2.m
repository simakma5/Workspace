% Experimental Data Analysis: Lecture 2
close all, clear all, clc

%% MATLAB EXAMPLE 1 %%
% generate data for normal distribution

data = randn(1,100); % generate random data for normal distribution

% make a figure
figure();
hist(data, 20);
ax = axis;
hold on

% plot a data
plot(data,0.3,'rx')

% set a figure properties
set(gca, 'Box', 'off')
set(gca, 'YLim', [0, 20])
set(gca, 'XLim', [-3, 3])
ylabel ('Frequency','FontSize',12, 'FontName','Georgia')
xlabel ('Value','FontSize',12, 'FontName','Georgia')

% compute basic data summary metrics
mn = mean(data);  % compute mean
sd = std(data);   % compute standard deviation
mabsdev = mad(data);    % compute mean absolute deviation
mabsdev = mad(data,1);    % compute median absolute deviation

% plot lines showing mean and +/- 1 std dev
h1 = plot([mn mn],      ax(3:4),'r-','LineWidth',2);
h2 = plot([mn-sd mn-sd],ax(3:4),'r-','LineWidth',2);
h3 = plot([mn+sd mn+sd],ax(3:4),'r-','LineWidth',2);
  % plot lines showing percentiles

%% generate data for non-normal distribution

data = (randn(1,100)).^2 + 1;  % generate random data for normal distribution

% make a figure
figure();
hist(data, 20);
ax = axis;
hold on

% plot a data
plot(data,0.3,'rx')

% set a figure properties
set(gca, 'Box', 'off')
set(gca, 'YLim', [0, 100])
set(gca, 'XLim', [0, 10])
ylabel ('Frequency','FontSize',12, 'FontName','Georgia')
xlabel ('Value','FontSize',12, 'FontName','Georgia')

% compute basic data summary metrics
med = median(data); % compute median
ptiles = prctile(data,[25 75]);  % compute 25th and 75th percentiles
% ptiles = iqr(data) % compute interquartile range, Q1-Q3, 25th-75th percentile

% plot lines showing mean and +/- 1 std dev
h1 = plot([med med],             ax(3:4),'r-','LineWidth',2);
h2 = plot([ptiles(1) ptiles(1)], ax(3:4),'r-','LineWidth',2);
h3 = plot([ptiles(2) ptiles(2)], ax(3:4),'r-','LineWidth',2);

%% MATLAB EXAMPLE 2 %%
% trimmed mean

data = (randn(1,100)).^2 + 1;  % generate random data for normal distribution

% make a figure
figure();
hist(data, 20);
ax = axis;
hold on

% plot a data
plot(data,0.3,'rx')

% set a figure properties
set(gca, 'Box', 'off')
set(gca, 'YLim', [0, 40])
set(gca, 'XLim', [0, 7])
ylabel ('Frequency','FontSize',12, 'FontName','Georgia')
xlabel ('Value','FontSize',12, 'FontName','Georgia')

% compute basic data summary metrics
med = median(data); % compute median
mn = mean(data);

% compute mean value while excluding the highets and lowest k data values,
% where k=n*(percent/100)/2
percent = 25;
trimmedmean = trimmean(data, percent); 

% plot lines showing mean and +/- 1 std dev
h1 = plot([med med], ax(3:4),'r-','LineWidth',2);
h2 = plot([mn mn], ax(3:4),'k-','LineWidth',2);
h3 = plot([trimmedmean trimmedmean], ax(3:4),'g-','LineWidth',2);

%% MATLAB EXAMPLE 3 %%
% generate probability density function
mn = 0;
sigma = 1;
X = [-4:0.1:4];
y = normpdf(X,mn,sigma);
plot(X,y,'r','LineWidth',2)
ylabel ('p(x)','FontSize',12, 'FontName','Georgia')
xlabel ('x','FontSize',12, 'FontName','Georgia')

%% MATLAB EXAMPLE 4 %%
% generate skewed distribution
% r = pearsrnd(mu,sigma,skew,kurt,m,n)
clc
data = pearsrnd(0,1,0,3,1000,1); % generate random data for normal distribution

% make a figure
figure();
hist(data, 20);
ax = axis;
hold on

% set a figure properties
set(gca, 'Box', 'off')
% set(gca, 'YLim', [0, 20])
% set(gca, 'XLim', [-3, 3])
ylabel ('Frequency','FontSize',12, 'FontName','Georgia')
xlabel ('Value','FontSize',12, 'FontName','Georgia')

skewness(data) % compute skewness
kurtosis(data) % compute curtosis

%% MATLAB EXAMPLE 5 %%
% graphical representation of data with normal distribution

subplot(121)
data = randn(1,100); % generate random data for normal distribution
[counts, bins] = hist(data, 20);
barh(bins,counts)% plot horizontal histogram

set(gca, 'XLim', [0, 20])
set(gca, 'YLim', [-3, 3])

subplot(122)
mn = mean(data);  % compute mean
sd = std(data);   % compute standard deviation
errorbar(10,mn,sd,'^','MarkerSize', 11,'Linestyle','--', 'linewidth', 1,'MarkerFaceColor','b','Color','k') % plot errorbar

set(gca, 'XLim', [0, 20])
set(gca, 'YLim', [-3, 3])

sem=std(data)/sqrt(length(data)); % compute standard error of the mean
ci = sem * 1.96; % compute 95% confidence interval

%% MATLAB EXAMPLE 6 %%
% graphical representation of data with non-normal distribution
data = (randn(1,100)).^2 + 1;  % generate random data for non-normal distribution

figure(), boxplot(data); % box plot

%% MATLAB EXAMPLE 7 %%
% bootstrapping
data = (randn(1,20)).^2 + 1;  % generate random data for non-normal distribution
hist(data) % plot a data
median(data)

m = bootstrp(1000,@median,data); % bootstrapping the median
figure(), hist(m); % plot a data
median(m)

%% MATLAB EXAMPLE 8 %%
% kernel density estimation
data = (randn(1,200)).^2 + 1;  % generate random data for normal distribution
X = [0:0.1:20]; % generate time vector
plot(data,0.05,'rx')
hold on
[f,x] = ksdensity(data,X); % estimation of probability distribution via kernel density estimation 
% f = f/sum(f); % normalize distribution to the area of 1
plot(X,f,'k','LineWidth',2) % plotting the resulting probability distribution 


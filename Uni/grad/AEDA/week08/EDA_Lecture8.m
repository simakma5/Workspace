% Experimental Data Analysis: Lecture 8
close all, clear all, clc

%% MATLAB EXAMPLE 1 %%
%  Fitting linearized regression model

% generate random data
x = 0:0.05:20;

% linearized model y = w1x^2 + w2x + w3
y = 2*x.^2 - 5*x + 10 + 20*randn(size(x));

% plot data
figure(),
plot(x, y,'x','MarkerSize',5,'LineWidth',2)
xlabel('x');
ylabel('y');

% regressor matrix
X = [x(:).^2 x(:) ones(length(x),1)];

% estimate the weights using ordinary least-squares
w = inv(X'*X)*X'*y(:)

% w = X\y(:); alternative approach

% estimate the model fit
modelfit = X*w;

% estimate the squared error
squarederror = sum((y(:)-modelfit).^2)

% plot the model
hold on;
plot(x,modelfit,'k-','LineWidth',2);

%%
% MATLAB EXAMPLE 1: Alternative approach using polyfit

close all

% fit the data x,y with the polynom of order 2
p = polyfit(x,y,2);

% compare w and p
w
p

% estimate the model fit
modelfit = polyval(p,x);

% plot the model
figure(),
plot(x, y,'x','MarkerSize',5,'LineWidth',2)
xlabel('x');
ylabel('y');
hold on;
plot(x,modelfit,'k-','LineWidth',2);

%% MATLAB EXAMPLE 2 %%
%  Fitting nonlinear regression model
close all, clear all, clc

% generate random data
x = 0:0.01:1;

% nonlinear model y = ax^n;
y = 4*x.^3 + 0.3*randn(size(x));

% plot data
figure(),
plot(x, y,'x','MarkerSize',5,'LineWidth',2)
xlabel('x');
ylabel('y');

% define optimization options
options = optimset('Display','iter','FunValCheck','on','MaxFunEvals',Inf,'MaxIter',Inf, 'TolFun',1e-6,'TolX',1e-6);

% define bounds for the parameters
%              a    n
paramslb = [-Inf    -Inf];  % lower bound
paramsub = [ Inf  Inf];  % upper bound

% define the initial seed
%            a     n
params0 = [  1     1];               

% estimate model parameters using nonlinear optimization
modelfun = @(par,data) par(1)*data.^par(2);

% estimate 'a' and 'n' parameters?
[params,resnorm,residual,exitflag,output] = lsqcurvefit(modelfun,params0,x,y,paramslb,paramsub,options);
params

% estimate the model fit
modelfit = modelfun(params,x);

% estimate the squared error
squarederror = sum((y(:)-modelfit(:)).^2)

% plot the model
hold on;
plot(x,modelfit,'k-','LineWidth',2);


%% MATLAB EXAMPLE 3 %%
%  Fitting nonlinear regression model

% stuttering (dysfluent words/second)
y=[3.983 8.044 10.136 4.547 6.047 1.713 7.150 8.424 0.377 1.817 2.296 1.506 2.525 5.392]';

% levodopa doses (mg)
x=[1210 1284 1402 503 1709 573 1101 928 683 336 577 238 619 305]';

% motor performance
x(:,2)=[4 15 6 29 27 11 33 23 36 16 22 28 27 30]';

% multiple linear regression
mdl = fitlm(x,y,'linear')

%% MATLAB EXAMPLE 4 %%
%  Comparing one-way ANOVA and linear regression
close all, clear all, clc

x = randi(3,100,1);     % Generate 100 random integer per 3 group
y = x + randn(100,1);   % Generate some response sample
x = categorical(x);  % treat the integers are categories/groups

% One-way ANOVA
[p,anovatab,stats] = anova1(y,x);

% Linear regression
mdl = fitlm(x,y);

% Show the group means from the ANOVA
ANOVAgroupMeans = stats.means

% Get the coeficients from the linear regression
LRcoefficients = [mdl.Coefficients.Estimate'] 

% Rescale the coeficients according the intercept
scaledLRcoefficients  = [LRcoefficients(1) LRcoefficients(1)+LRcoefficients(2:3)]

% Compare the two results
abs(max(scaledLRcoefficients - ANOVAgroupMeans)) 

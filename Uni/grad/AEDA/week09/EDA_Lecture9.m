% Experimental Data Analysis: Lecture 9
close all, clear all, clc

%% MATLAB EXAMPLE 1 %%
%  Logistic regression

load data;
X = data; % features

y(1:40) = 0; % the first half of data belonging to group 1
y(41:80) = 1; % the second half of data belonging to group 2

% visualize the data
figure(); hold on;
h1 = plot(X(y==0,1),X(y==0,2),'b.','MarkerSize',20);
h2 = plot(X(y==1,1),X(y==1,2)','ko','MarkerSize',6);
legend([h1 h2],{'Huntington' 'Healthy'},'Location','NorthEast');
xlabel('Articulation rate');
ylabel('Vowel articulation');

% perform a logistic regression
modelparams = glmfit(X,y','binomial','link','logit');

% compute the class assignment (probabilities) of the trained model for each data point
modelfitLR = glmval(modelparams,X,'logit') >= 0.5;

% estimate final classification accuracy
LRcorrect = sum(modelfitLR==y')/length(y)*100

% create a grid matrix to visualise the model
ax = axis; % automatically compute handles
xvals = linspace(ax(1),ax(2),1000); % range of 1. feature values
yvals = linspace(ax(3),ax(4),1000); % range of 2. feature values
[xx,yy] = meshgrid(xvals,yvals); % create a grid for image
gridX = [xx(:) yy(:)]; % transform grid into N x number of features vector suitable for classification

% perform a logistic regression with all grid points
outputimageLR = glmval(modelparams,gridX,'logit');

% transform grid back to image
outputimageLR = reshape(outputimageLR,[length(yvals) length(xvals)]);

% draw a decision boundary line
[~,h3] = contour(xvals,yvals,outputimageLR,[.5 .5]); % decision at the point 0.5
set(h3,'LineWidth',2,'LineColor',[0 0 0]);  
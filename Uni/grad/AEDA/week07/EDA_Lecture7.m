% Experimental Data Analysis: Lecture 7
close all, clear all, clc

%% MATLAB EXAMPLE 1 %%
% Two-way ANOVA

% load data  
data = [
1	1	5
1	1	7
1	1	8
1	2	5
1	2	5
1	2	4
2	1	8
2	1	17
2	1	11
2	2	10
2	2	15
2	2	10
3	1	29
3	1	22
3	1	18
3	2	18
3	2	13
3	2	11
];

group = data(:,1); % 1st Factor: Group
sex = data(:,2); % 2nf Factor: Sex
sample = data(:,3); % Individual scores

% calculate two-way ANOVA
[p table stats terms] = anovan(sample, {group, sex},'varnames',{'group', 'sex'},'model','interaction');

% calculate post-hoc differences
[c,m,h,nms] = multcompare(stats,'dimension',1,'ctype','bonferroni','alpha',0.05);
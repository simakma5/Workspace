% Experimental Data Analysis: Lecture 11
clear all, close all, clc

%% MATLAB EXAMPLE 1 %%
%  Clustering
load data;

%% K-means
%     F1a F2a   F1  F2i   F1u F2u
mu = [700 1500; 350 2300; 350 800]; % define starting positions
options = statset('Display','final'); % define display options

% perform k-means
[objy, ctrs] = kmeans(data,3,'Options',options,'Start',mu); 

% assign found centroids to formant frequencies
F1a = ctrs(1,1);
F2a = ctrs(1,2);
F1i = ctrs(2,1);
F2i = ctrs(2,2);
F1u = ctrs(3,1);
F2u = ctrs(3,2);

% plot vowel space triangle
subplot(121),
plot(data(:,1),data(:,2),'kx');
hold on
plot(F1a,F2a,'rx','markersize',14,'LineWidth',4)
plot(F1i,F2i,'gx','markersize',14,'LineWidth',4)
plot(F1u,F2u,'bx','markersize',14,'LineWidth',4)
set(gca,'box','off')

F = [F1a F2a; F1i F2i; F1u F2u];
DT = delaunayTriangulation(F(:,1),F(:,2));
[k, v] = convexHull(DT);
plot(DT.Points(k,1),DT.Points(k,2),'c','LineWidth',3);

%% EM algorithm

%       F1a F2a   F1  F2i   F1u F2u
S.mu = [700 1500; 350 2300; 350 800]; % define starting positions
S.Sigma = cat(3,[100 0; 0 100],[100 0;0 100],[100 0;0 100]); % define variances sigma
S.PComponents = ones(1,3)/3; % define Fí
options = statset('Display','final'); % define display options

% run EM algorithm
objx = fitgmdist(data,3,'Options',options,'Start',S);
% objx = gmdistribution.fit(fX,5,'Options',options,'Start','randsample','Replicates',10);

% assign found centroids to formant frequencies
F1a = objx.mu(1,1);
F2a = objx.mu(1,2);
F1i = objx.mu(2,1);
F2i = objx.mu(2,2);
F1u = objx.mu(3,1);
F2u = objx.mu(3,2);

% plot vowel space triangle
subplot(122),
plot(data(:,1),data(:,2),'kx');
hold on
plot(F1a,F2a,'rx','markersize',14,'LineWidth',3)
plot(F1i,F2i,'gx','markersize',14,'LineWidth',4)
plot(F1u,F2u,'bx','markersize',14,'LineWidth',4)
set(gca,'box','off')

F = [F1a F2a; F1i F2i; F1u F2u];
DT = delaunayTriangulation(F(:,1),F(:,2));
[k, v] = convexHull(DT);
plot(DT.Points(k,1),DT.Points(k,2),'c','LineWidth',4);
